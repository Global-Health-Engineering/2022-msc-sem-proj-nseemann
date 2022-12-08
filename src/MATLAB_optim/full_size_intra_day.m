%%%ETH ZURICH 
%%%Course: SP
close all
%clear

% Change to current running directory:
% cd 'D:\Semester_project_git\2022-msc-sem-proj-nseemann\src\MATLAB_optim'

%% Parameters
% params.optiParams.period_op_cost = 4.5;
% params.optiParams.truck_cost = 0;
% params.optiParams.km_cost = 0.5;
% params.optiParams.skip_add_cost = 0;
% params.optiParams.max_add_bins = 13;
% params.optiParams.period_t_max = 4;
% 
% params.optiParams.numV = 3;
% params.preParams.speed_avg = 30;
% params.preParams.dump_ind = 54;
% params.preParams.compost_ind = 55;
% params.preParams.depot_ind = 56;
% params.preParams.T = 7;
% params.preParams.P = 2;
% params.preParams.underfull_threshold = 0.6;
% params.preParams.set_add_bins = 2;
% params.preParams.per_week_consider = 3;
% params.preParams.unified_salary = 1;
%%

% % Cost of operating one truck for 1 day (assume labour costs)
% period_op_cost = 4.5;
% 
% % Average speed 30km/h
% speed_avg = 30; 
% 
% % Capital cost of truck / maximum number of trucks
% truck_cost = 0; % 20000 according to Liz;
% 
% %max_truck = 4;
% numV = 3;
% % Distance cost
% km_cost = 0.5;
% 
% % Capital cost of buying one skip
% skip_add_cost = 0;% 3000 according to Liz;
% max_add_bins = 10;
% 
% % Indices other than skips
% dump_ind = 54; % Mzedi dump
% compost_ind = 55; % Limbe composting
% depot_ind = 56; % Storage city center


% Time days
T = params.preParams.T;
% Number of collections periods per day
P = params.preParams.P;

% % Fulness threshold under which scenarios are not considered
% underfull_threshold = 0.6;
% 
% % maximum number of additional bins at each bin
% set_add_bins = 2;
% 
% % additional bins only considered for more or equal than per_week_consider services
% % under current service
% per_week_consider = 4;
% 
% period_t_max = 4; % Maximum number of hours per period
% 
% unified_salary = 1; % If the wage is constant for all the crews no matter when the crew is dispatched


%% Import distances
%Distance matrix
fileID = fopen('../../data/interm_data/dist_matrix_ID_filtered.csv');
textscan(fileID,'%s',3,'Delimiter',',');
C_dist = textscan(fileID,'%d %d %f','Delimiter',',');
fclose(fileID);


%Convert into a matrix
dist_mat = zeros(length(max(max(C_dist{1},max(C_dist{2})))));
for i = 1:length(C_dist{1})
    dist_mat(C_dist{1}(i),C_dist{2}(i)) = C_dist{3}(i);
end

%Travel time matrix (assuming constant speed)
t_mat = dist_mat/params.preParams.speed_avg;

%% Filtering distances and isolating dump, composting and storage facility
indices_facilities = sort([params.preParams.dump_ind,params.preParams.compost_ind,params.preParams.depot_ind]); % Indices of facilities sorted

% Indices of skips
indices_skips = [1:indices_facilities(1)-1,indices_facilities(1)+1:indices_facilities(2)-1,indices_facilities(2)+1:indices_facilities(3)-1,indices_facilities(3)+1:max(max(C_dist{1},max(C_dist{2})))];

% Number of skips
numBins = length(indices_skips);

%% Import filling_rates
    %Filling rates should have dimension 1 x numBins or numBins x T*P for
    %variable
fileID = fopen('../../data/interm_data/skips_indexed_filling_estimates.csv');
C_rates_headings = textscan(fileID,'%s',4,'Delimiter',',');
C_rates = textscan(fileID,'%s %d %f %f','Delimiter',',');
fclose(fileID);
% Random but constant rates
% load('filling_rates_13_10_22.mat')
filling_rates = zeros(numBins,T*P);
for i = 1:numBins
    filling_rates(indices_skips(i),:) = repmat(C_rates{4}(find(C_rates{2} == indices_skips(i))),1,T*P);
end
filling_rates(isnan(filling_rates)) = params.preParams.nanreplace;
%filling_rates = rand(1,numBins)/2;
%histogram(filling_rates)
%ylabel('number of skips')
%xlabel('filling rate [/day]')


%% Load scenarios and relevant variables
% Load scenarios and largest gap in scenario
[scens,mult]=create_scens(); %Week multiplier, for extra-weekly collections
divis = 1./mult; %Divisor for costs used in objective function

%% Narrow relevant scenarios
scen_cells = {}; %Cell for each skip. {{scen numbers}, {number of additional skips}, {cost multiplier}} 
%if a scenarios is suitable, add the scenario vector to skip scenario cell
for i = 1:numBins
    scen_cells{i,1} = [];
    scen_cells{i,2} = [];
    scen_cells{i,3} = [];
    
    % Place the rate adjacent to itself to get entire cycle filling
    current_rates_doubled = [filling_rates(i,:) filling_rates(i,:)];
    
    % For each scenario
    for s = 1:size(scens,1)
        
        % Places the scenarios adjacent to itself to get the entire cycle
        current_scen_doubled = [scens(s,:) scens(s,:)]; 
        
        previous = 1; % Previous index 
        exceed = 0; % Flag if the sum exceeds
        
        fullness_levels = [];
        % For all indices in the scenario where the bin is emptied
        for d = find(current_scen_doubled) 
            % Accumulated waste between previous index and current one,
            % times filling multiplier for extra weekly schedules
            sum_period = sum(current_rates_doubled(previous:d)*mult(s))/P;
            if sum_period >= 1 % If skip overfull
                exceed = 1; % Raise flag
                break
            end
            fullness_levels = [fullness_levels sum_period];
            previous = d + 1; % Move previous cursor to current + 1
        end 
        % If the exceed flag is not raised and if not (intra weekly undefull) i.e.
        % underfull extra_weekly are allowed
        if exceed ~= 1 && ~(sum(scens(s,:)) > 1 && max(fullness_levels) + min(filling_rates(i,:)) < params.preParams.underfull_threshold)
            scen_cells{i,1} = [scen_cells{i,1}; s]; % Add scen number to first cell column
            scen_cells{i,2} = [scen_cells{i,2}; 0]; % 0 additional skips for now
            scen_cells{i,3} = [scen_cells{i,3}; divis(s)]; % cost divisor
        end
    end
    
    % if extra weekly, only keep the biggest one and delete all other
    % scenarios
    max_extra_period = 1./min(scen_cells{i,3});
    if max_extra_period > 1
        scen_cells{i,1} = scen_cells{i,1}(scen_cells{i,3} == min(scen_cells{i,3}));
        scen_cells{i,2} = scen_cells{i,2}(scen_cells{i,3} == min(scen_cells{i,3}));
        scen_cells{i,3} = scen_cells{i,3}(scen_cells{i,3} == min(scen_cells{i,3}));
    end

    % Adding scenarios for added bins up to maximum number of allowable
    % additional bins per bin location
    for add_bins = 1:params.preParams.set_add_bins   
        %Scens with more than per_week_consider in a week or if possible
        %schedules is empty (extra weekly underfull have the largest extra
        %weekly scenario assigned above
        if isempty(scen_cells{i,1}) || min(sum(scens(scen_cells{i,1},:),2).*scen_cells{i,3}) >= params.preParams.per_week_consider 
            current_rates_doubled = [filling_rates(i,:)./(1+add_bins) filling_rates(i,:)./(1+add_bins)]; %doubled for extra weekly
            for s = 1:size(scens,1)
                current_scen_doubled = [scens(s,:) scens(s,:)]; %doubled for extra weekly
                previous = 1;
                exceed = 0;
                fullness_levels = [];
                for d = find(current_scen_doubled)
                    sum_period = sum(current_rates_doubled(previous:d)*mult(s))/P;
                    if sum_period >= 1
                        exceed = 1;
                        break
                    end
                    fullness_levels = [fullness_levels sum_period];
                    previous = d + 1;
                end 
                if exceed ~= 1 && ~(sum(scens(s,:)) > 1 && max(fullness_levels) < params.preParams.underfull_threshold)
                    scen_cells{i,1} = [scen_cells{i,1}; s]; %scens(s,:) <-- old entire scen
                    scen_cells{i,2} = [scen_cells{i,2}; add_bins];
                    scen_cells{i,3} = [scen_cells{i,3}; divis(s)];
                end
            end
        end
    end 
end

possible_extra = 0;
for i = 1:length(scen_cells)
    possible_extra = possible_extra + max(scen_cells{i,2});
end

scen_cells_extra = {};
skip_nums_extra_scens = [];
scen_cells_indices = {};
for l = 1:size(scen_cells,1) % For all scens 
    for i = 1:max(scen_cells{l,2}) %For each additional bin number
        for j = 1:i %Run as many times as there are bins
            scen_cells_extra = [scen_cells_extra scen_cells{l,1}(scen_cells{l,2}==i)]; %Possible scenarios for extra skips
            skip_nums_extra_scens = [skip_nums_extra_scens l]; %Bin corresponding to extra bin
            scen_cells_indices = [scen_cells_indices find(scen_cells{l,2} == i)]; %Position in yis
        end
    end
end

t_mat_dump_extra = t_mat(params.preParams.dump_ind,skip_nums_extra_scens);
t_mat_depot_extra = t_mat(params.preParams.depot_ind,skip_nums_extra_scens);
%% Decision variables and constraints

% Scenario selection for regular skips
yis = {};
for l = 1:size(scen_cells,1)
    yis{l} = binvar(size(scen_cells{l,1},1), 1,'full'); %if l is assigned scenario in subset
end


%%
%Flow 

xit = binvar(T*P,numBins,'full'); %if is operating on day t at period p
xit_extra = binvar(T*P,size(scen_cells_extra,2));
%% Decision variables

fit = binvar(T*P, numBins,'full'); %if first on day t at p
fit_extra = binvar(T*P,size(scen_cells_extra,2),'full');

numV_D = intvar(T*P,1,'full');


if params.preParams.unified_salary == 0
    numV_crews = intvar(T*P,1,'full'); % Number of vehicles in operation on day
else 
    numV_crews = [repmat(intvar(1,1),T*P,1)];
end

%numV = intvar(1); % Total number of vehicles


%% Constraints
% Equality constraints

assign_1 = []; %Only one scenario within allowed scenario per skip
for l = 1:size(scen_cells,1)
    assign_1 = [assign_1 sum(yis{l})==1];
end

yis_extra = {};
assign_1_extra = [];
for l = 1:size(scen_cells_extra,2) % For each extra skip
    yis_extra{l} = binvar(size(scen_cells_extra{l},1), 1,'full');
    assign_1_extra  = [assign_1_extra sum(yis_extra{l}) == 1];
end

% NORMAL CONSTRAINT HERE FOR SCENARIO SELECTION EXTRA

day_flow = []; %Skips emptied only on periods assigned by scenario
for l = 1:size(scen_cells,1)
    day_flow = [day_flow xit(:,l)'-yis{l}'*scens(scen_cells{l,1},:)==0];
end

disp('Day flow extra starting');
day_flow_extra = []; %Skips emptied only on periods assigned by scenario
for l = 1:size(scen_cells_extra,2)
    apply_mult = sum(yis{skip_nums_extra_scens(l)}(scen_cells_indices{l}));
    for i = 1:T*P
        day_flow_extra = [day_flow_extra xit_extra(i,l)'-apply_mult*yis_extra{l}'*scens(scen_cells_extra{l},i)==0];
    end
    toc;
end
disp('Day flow extra done');


% Inequality constraints
crews = numV_crews >= numV_D;
% Ensure period time is not exceeded %MUST ADD FIT
%   period_time= xit_extra*t_mat(params.preParams.dump_ind,skip_nums_extra_scens)' + (t_mat(skip_nums_extra_scens,params.preParams.dump_ind)'*xit_extra')' +  xit*t_mat(params.preParams.dump_ind,indices_skips)' + (t_mat(indices_skips,params.preParams.dump_ind)'*xit')' + fit*((-1*t_mat(params.preParams.dump_ind,indices_skips) + t_mat(params.preParams.depot_ind,indices_skips)))' <= repmat(params.optiParams.period_t_max,T*P,1).*numV_D;

period_time = [xit_extra xit]*t_mat(params.preParams.dump_ind,[skip_nums_extra_scens indices_skips])'...
    + (t_mat([skip_nums_extra_scens indices_skips],params.preParams.dump_ind)'*[xit_extra xit]')'...
    + [fit_extra fit]*((-1*t_mat(params.preParams.dump_ind,[skip_nums_extra_scens indices_skips]) + t_mat(params.preParams.depot_ind,[skip_nums_extra_scens indices_skips])))'...
    + numV_D*t_mat(params.preParams.dump_ind,params.preParams.depot_ind) <= repmat(params.optiParams.period_t_max,T*P,1).*numV_D;

% period_time = [xit_extra xit]*t_mat(params.preParams.dump_ind,[skip_nums_extra_scens indices_skips])'...
%     + (t_mat([skip_nums_extra_scens indices_skips],params.preParams.dump_ind)'*[xit_extra xit]')'...
%     + [fit_extra fit]*((-1*t_mat(params.preParams.dump_ind,[skip_nums_extra_scens indices_skips]) + t_mat(params.preParams.depot_ind,[skip_nums_extra_scens indices_skips])))'...
%     + sum([fit_extra fit],2)*t_mat(params.preParams.dump_ind,params.preParams.depot_ind) <= repmat(params.optiParams.period_t_max,T*P,1).*numV_D;

%Fit applying to fit_extra
first = sum(fit, 2) + sum(fit_extra,2) == numV_D; %Makes it so numV_D skip is first in loop on each period
fit_unity = fit <= xit; %Prevents a first in loop when a skip is not operating on that period
fit_extra_unity = fit_extra <= xit_extra;

%Number of trucks under maximum 
%numV_const = 1<= numV <= max_truck;

%Number of trucks per day under overall number of trucks
numV_D_numV = zeros(T*P,1) <= numV_D <= repmat(params.optiParams.numV,T*P,1);

total_added_skips = 0;
for l = 1:size(scen_cells,1)
   total_added_skips = total_added_skips + yis{l}'*scen_cells{l,2};
end

add_skip_max = total_added_skips <= params.optiParams.max_add_bins;
% Overall constraints
Constraints = [crews assign_1 fit_extra_unity day_flow_extra assign_1_extra  numV_D_numV day_flow  first fit_unity period_time  add_skip_max];
%fill_min
%period_op
%numV_const
%% Objective function

% Divisors by skip
divis_vect = [];
for l = 1:size(scen_cells,1)
   divis_vect = [divis_vect yis{l}'*scen_cells{l,3}];
end

add_bins_vect = [];
for l = 1:size(scen_cells,1)
   add_bins_vect = [add_bins_vect yis{l}'*scen_cells{l,2}];
end


for l = 1:size(scen_cells,1)
    if(sum(scen_cells{l,3}) == length(scen_cells{l,3}))
        xit_adj(:,l) = xit(:,l);
        fit_adj(:,l) = fit(:,l);
    else
        xit_adj(:,l) = xit(:,l).*scen_cells{l,3}(1);
        fit_adj(:,l) = fit(:,l); %Don't divide fit, to not underestimate costs
    end
end
%%
% Distance travelled on the rounds
% distance_diff_old = params.optiParams.km_cost*( xit_extra*dist_mat(params.preParams.dump_ind,skip_nums_extra_scens)' + (dist_mat(skip_nums_extra_scens,params.preParams.dump_ind)'*xit_extra')' + xit_adj*(dist_mat(params.preParams.dump_ind,indices_skips))' + (dist_mat(indices_skips,params.preParams.dump_ind)'*xit_adj')' + fit_adj*((-1*dist_mat(params.preParams.dump_ind,indices_skips) + dist_mat(params.preParams.depot_ind,indices_skips)))');
distance_cost = params.optiParams.km_cost*(...
    [xit_extra xit_adj]*dist_mat(params.preParams.dump_ind,[skip_nums_extra_scens indices_skips])'...
    + (dist_mat([skip_nums_extra_scens indices_skips],params.preParams.dump_ind)'*[xit_extra xit_adj]')'...
    + [fit_extra fit_adj]*((-1*dist_mat(params.preParams.dump_ind,[skip_nums_extra_scens indices_skips]) + dist_mat(params.preParams.depot_ind,[skip_nums_extra_scens indices_skips])))'...
    + numV_D*dist_mat(params.preParams.dump_ind,params.preParams.depot_ind));
% Daily operation cost
%capital_cost =  params.optiParams.skip_add_cost*sum(add_bins_vect) + params.optiParams.truck_cost*params.optiParams.numV;
total_op_cost = params.optiParams.period_op_cost.*numV_crews; %sum(ot)
Objective = sum(distance_cost) + sum(total_op_cost); %+ capital_cost;

%% Set options for YALMIP and solver - Solve the problem
% CUTSDP or gurobi
options =   params.optiOptions;%sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1,'gurobi.MIPGap',0.05);
disp('optimize being called');
sol =       optimize(Constraints,Objective,options);
disp('optimized finished');
toc
%% Analyze error flags & Get the solution
if sol.problem == 0  
    disp('Finished running');

    
else
    disp('Hmm, something went wrong!');
    sol.info;
    yalmiperror(sol.problem)
end

toc;

%% Solution processing

output.process.scen_cells = scen_cells;
output.process.possible_extra = possible_extra;
output.process.scen_cells_extra = scen_cells_extra;
output.process.scen_cells_indices = scen_cells_indices;
output.process.skip_nums_extra_scens = skip_nums_extra_scens;
output.process.scens = scens;
output.process.C_rates = C_rates;
output.process.filling_rates = filling_rates;
output.process.dist_mat = dist_mat;
output.process.t_mat = t_mat;
output.process.indices_facilities = indices_facilities;
output.process.indices_skips = indices_skips;

if ~strcmp(sol.solveroutput.result.status,'INF_OR_UNBD')
    output.optiVars.xit = value(xit);
    output.optiVars.xit_adj = value(xit_adj);
    output.optiVars.xit_extra = value(xit_extra);
    output.optiVars.fit = value(fit);
    output.optiVars.fit_extra = value(fit_extra);

    for i = 1:length(yis)
        output.optiVars.yis{i} = value(yis{i});
    end

    for i = 1:length(yis_extra)
        output.optiVars.yis_extra{i} = value(yis_extra{i});
    end
    output.optiVars.numV_D = value(numV_D);
    output.optiVars.numV_crews = value(numV_crews);

    output.results.distance_cost = value(distance_cost);
    output.results.total_op_cost = value(total_op_cost);
    output.results.add_bins_vect = value(add_bins_vect);
    output.results.Objective = value(Objective);



    output.solver_output = sol;

    output.bound_gap = (sol.solveroutput.result.objval-sol.solveroutput.result.poolobjbound)/sol.solveroutput.result.objval;
    output.feasible = 1;
else
    output.feasible = 0;
end
%% Builds scenarios
%1x14 double period one weekly with week multiplier/cost divisor
function [scens,mult] = create_scens()
    reserved_periods = [13 14]; %Sundays
    P = 2;
    T=7;
    % Every period (twice a day) except Sunday
    scens = [1 1 1 1 1 1 1 1 1 1 1 1 0 0];
    mult = [1];
    %once every 4,3,2,1 weeks:
    num_new_scens = 0;
    for j = 1:4
        for i = 1:T*P
         if ~ismember(i,reserved_periods)
          scens = [scens;repmat([0],1,i-1) 1 repmat([0], 1,(T*P)-i)];
          num_new_scens = num_new_scens + 1;
         end
        end
    end
    for j = 1:4
        for i = 1:num_new_scens/4
            mult = [mult j];
        end
    end
    
    %disp(scens);
    
    scens_intra_week = [];
    basis_D = [1 0 0 0 0 0 0;0 1 0 0 0 0 0; 0 0 1 0 0 0 0;0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0];  
    basis_D = [basis_D;1 0 0 1 0 0 0;0 1 0 0 0 1 0; 1 0 0 0 1 0 0;0 0 1 0 0 1 0; 0 1 0 0 1 0 0 ]; %2/week
    basis_D = [basis_D;1 0 1 0 1 0 0; 0 1 0 1 0 1 0; 1 0 0 1 0 1 0; 1 0 1 0 0 1 0]; %3/week
    basis_D = [basis_D; 1 0 1 0 1 1 0; 1 1 0 1 0 1 0; 1 0 1 1 0 1 0]; %4/week
    basis_D = [basis_D; 1 0 1 1 1 1 0; 1 1 0 1 1 1 0; 1 1 1 0 1 1 0; 1 1 1 1 0 1 0]; %5/week
    basis_D = [basis_D; 1 1 1 1 1 1 0;]; %6/week
    scen_built = [];
    for i =  1:size(basis_D,1)
        find_basis = find(basis_D(i,:));
        scen_doubled = [];
        for j = 1:size(basis_D(i,:),2)
            scen_doubled = [scen_doubled basis_D(i,j) basis_D(i,j)]; 
        end
        cases_bin = dec2bin(0:2^(size(find_basis,2))-1) - '0';
        for j = 1:size(cases_bin,1)
            scen_built = scen_doubled;
            for k = 1:size(cases_bin,2)
                scen_built((find_basis(k)*2)-1:(find_basis(k)*2)) = [cases_bin(j,k) ~cases_bin(j,k)];
            end
            scens_intra_week = [scens_intra_week; scen_built];
        end
    end
    rem_ind = [];
    % Filter scens so that maxgap - mingap <= 2
    for j = 1:size(scens_intra_week,1)
        diff_indices = diff(find([scens_intra_week(j,:) scens_intra_week(j,:)]));
        if max(diff_indices) - min(diff_indices) > 2
            rem_ind = [rem_ind j];
        end
    end
    scens_intra_week(rem_ind,:) = [];
    
    mult = [mult ones(1,size(scens_intra_week,1))];
    scens = [scens;scens_intra_week];
end 