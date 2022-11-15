%%%ETH ZURICH 
%%%Course: SP
close all
clear; 
tic;
disp('timer started');

% Change to current running directory:
% cd 'D:\Semester_project_git\2022-msc-sem-proj-nseemann\src\MATLAB_optim'

%% Parameters

% Cost of operating one truck for 1 day (assume labour costs)
period_op_cost = 4.5;

% Average speed 30km/h
speed_avg = 30; 

% Capital cost of truck / maximum number of trucks
truck_cost = 0; % 20000 according to Liz;
max_truck = 2;

% Distance cost
km_cost = 0.5;

% Capital cost of buying one skip
skip_add_cost = 0;% 3000 according to Liz;
max_add_bins = 40;

% Indices other than skips
dump_ind = 54; % Mzedi dump
compost_ind = 55; % Limbe composting
depot_ind = 56; % Storage city center


% Time days
T = 7;
% Number of collections periods per day
P = 2;

% maximum number of additional bins at each bin
set_add_bins = 3;

% additional bins only considered for more or equal than per_week_consider services
% under current service
per_week_consider = 5;

period_t_max = 4; % Maximum number of hours per period

%% Import distances
%Distance matrix
fileID = fopen('../../data/interm_data/dist_matrix_ID_filtered.csv');
textscan(fileID,'%s',3,'Delimiter',',');
C = textscan(fileID,'%d %d %f','Delimiter',',');
fclose(fileID);

%Convert into a matrix
dist_mat = zeros(length(max(max(C{1},max(C{2})))));
for i = 1:length(C{1})
    dist_mat(C{1}(i),C{2}(i)) = C{3}(i);
end

%Travel time matrix (assuming constant speed)
t_mat = dist_mat/speed_avg;

%% Filtering distances and isolating dump, composting and storage facility
indices_facilities = sort([dump_ind,compost_ind,depot_ind]); % Indices of facilities sorted

% Indices of skips
indices_skips = [1:indices_facilities(1)-1,indices_facilities(1)+1:indices_facilities(2)-1,indices_facilities(2)+1:indices_facilities(3)-1,indices_facilities(3)+1:max(max(C{1},max(C{2})))];

% Number of skips
numBins = length(indices_skips);

%% Import filling_rates
    %Filling rates should have dimension 1 x numBins or numBins x T*P for
    %variable
 
% Random but constant rates
load('filling_rates_13_10_22.mat')

%filling_rates = rand(1,numBins)/2;
%histogram(filling_rates)
%ylabel('number of skips')
%xlabel('filling rate [/day]')

% Makes uniform filling rates over entire period
filling_rates = repmat(filling_rates',1,T*P);

% Forcing some filling rates to extreme values
filling_rates(1,:) = [repmat([0.86],1, 14)];
filling_rates(5,:) = [repmat([0.7],1, 14)];
filling_rates(8,:) = [repmat([1.2],1, 14)];
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
        
        % For all indices in the scenario where the bin is emptied
        for d = find(current_scen_doubled) 
            
            % Accumulated waste between previous index and current one,
            % times filling multiplier for extra weekly schedules
            sum_period = sum(current_rates_doubled(previous:d)*mult(s))/P;
            if sum_period >= 1 % If skip overfull
                exceed = 1; % Raise flag 
                break
            end
            previous = d + 1; % Move previous cursor to current + 1
        end 
        if exceed ~= 1 % If exceed flag not raised
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
    for add_bins = 1:set_add_bins   
        %Scens with more than per_week_consider in a week or if possible
        %schedules is empty
        if isempty(scen_cells{i,1}) || min(sum(scens(scen_cells{i,1},:),2).*scen_cells{i,3}) >= per_week_consider 
            current_rates_doubled = [filling_rates(i,:)./(1+add_bins) filling_rates(i,:)./(1+add_bins)]; %doubled for extra weekly
            for s = 1:size(scens,1)
                current_scen_doubled = [scens(s,:) scens(s,:)]; %doubled for extra weekly
                previous = 1;
                exceed = 0;
                for d = find(current_scen_doubled)
                    sum_period = sum(current_rates_doubled(previous:d)*mult(s))/P;
                    if sum_period >= 1
                        exceed = 1;
                        break
                    end
                    previous = d + 1;
                end 
                if exceed ~= 1
                    scen_cells{i,1} = [scen_cells{i,1}; s]; %scens(s,:) <-- old entire scen
                    scen_cells{i,2} = [scen_cells{i,2}; add_bins];
                    scen_cells{i,3} = [scen_cells{i,3}; divis(s)];
                end
            end
        end
    end 
end

%% Decision variables and constraints
yis = {};
for l = 1:size(scen_cells,1)
    yis{l} = binvar(size(scen_cells{l,1},1), 1,'full'); %if l is assigned scenario in subset
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

t_mat_dump_extra = t_mat(dump_ind,skip_nums_extra_scens);
t_mat_depot_extra = t_mat(depot_ind,skip_nums_extra_scens);


yis_extra = {};
assign_1_extra = [];
for l = 1:size(scen_cells_extra,2) % For each extra skip
    yis_extra{l} = binvar(size(scen_cells_extra{l},1), 1,'full');
    assign_1_extra  = [assign_1_extra sum(yis_extra{l}) == 1];
end

%%
%Flow 
disp('Day flow extra starting');
xit_extra = binvar(T*P,size(scen_cells_extra,2));
day_flow_extra = []; %Skips emptied only on periods assigned by scenario
for l = 1:size(scen_cells_extra,2)
    apply_mult = sum(yis{skip_nums_extra_scens(l)}(scen_cells_indices{l}));
    for i = 1:T*P
        day_flow_extra = [day_flow_extra xit_extra(i,l)'-apply_mult*yis_extra{l}'*scens(scen_cells_extra{l},i)==0];
    end
    toc;
end
disp('Day flow extra done');
%% Decision variables
%matching decision var
numExtraBinScens = size(xit_extra,2);

xit = binvar(T*P,numBins,'full'); %if is operating on day t at period p
fit = binvar(T*P, numBins,'full'); %if first on day t at p
ot = binvar(T*P,1,'full'); % if operating on day t at p
numV_D = intvar(T*P,1,'full'); % Number of vehicles in operation on day

numV = intvar(1); % Total number of vehicles

fit_extra = binvar(T*P,numExtraBinScens,'full');
%% Constraints
% Equality constraints

assign_1 = []; %Only one scenario within allowed scenario per skip
for l = 1:size(scen_cells,1)
    assign_1 = [assign_1 sum(yis{l})==1];
end
% NORMAL CONSTRAINT HERE FOR SCENARIO SELECTION EXTRA



day_flow = []; %Skips emptied only on periods assigned by scenario
for l = 1:size(scen_cells,1)
    day_flow = [day_flow xit(:,l)'-yis{l}'*scens(scen_cells{l,1},:)==0];
end

% Inequality constraints

% Ensure period time is not exceeded %MUST ADD FIT
%   period_time= xit_extra*t_mat(dump_ind,skip_nums_extra_scens)' + (t_mat(skip_nums_extra_scens,dump_ind)'*xit_extra')' +  xit*t_mat(dump_ind,indices_skips)' + (t_mat(indices_skips,dump_ind)'*xit')' + fit*((-1*t_mat(dump_ind,indices_skips) + t_mat(depot_ind,indices_skips)))' <= repmat(period_t_max,T*P,1).*numV_D;

period_time = [xit_extra xit]*t_mat(dump_ind,[skip_nums_extra_scens indices_skips])'...
    + (t_mat([skip_nums_extra_scens indices_skips],dump_ind)'*[xit_extra xit]')'...
    + [fit_extra fit]*((-1*t_mat(dump_ind,[skip_nums_extra_scens indices_skips]) + t_mat(depot_ind,[skip_nums_extra_scens indices_skips])))'...
    + numV_D*t_mat(dump_ind,depot_ind) <= repmat(period_t_max,T*P,1).*numV_D;
%Fit applying to fit_extra
first = sum(fit, 2) + sum(fit_extra,2) == numV_D; %Makes it so numV_D skip is first in loop on each period
fit_unity = fit <= xit; %Prevents a first in loop when a skip is not operating on that period
fit_extra_unity = fit_extra <= xit_extra;

%Number of trucks under maximum 
numV_const = 1<= numV <= max_truck;

%Number of trucks per day under overall number of trucks
numV_D_numV = zeros(T*P,1) <= numV_D <= repmat(numV,T*P,1);

total_added_skips = 0;
for l = 1:size(scen_cells,1)
   total_added_skips = total_added_skips + yis{l}'*scen_cells{l,2};
end

add_skip_max = total_added_skips <= max_add_bins;
% Overall constraints
Constraints = [assign_1 fit_extra_unity day_flow_extra assign_1_extra  numV_D_numV day_flow  first fit_unity period_time numV_const add_skip_max];
%fill_min
%period_op
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
        fit_adj(:,l) = fit(:,l).*scen_cells{l,3}(1);
    end
end
%%
% Distance travelled on the rounds
% distance_diff_old = km_cost*( xit_extra*dist_mat(dump_ind,skip_nums_extra_scens)' + (dist_mat(skip_nums_extra_scens,dump_ind)'*xit_extra')' + xit_adj*(dist_mat(dump_ind,indices_skips))' + (dist_mat(indices_skips,dump_ind)'*xit_adj')' + fit_adj*((-1*dist_mat(dump_ind,indices_skips) + dist_mat(depot_ind,indices_skips)))');
distance_diff = km_cost*(...
    [xit_extra xit_adj]*dist_mat(dump_ind,[skip_nums_extra_scens indices_skips])'...
    + (dist_mat([skip_nums_extra_scens indices_skips],dump_ind)'*[xit_extra xit_adj]')'...
    + [fit_extra fit_adj]*((-1*dist_mat(dump_ind,[skip_nums_extra_scens indices_skips]) + dist_mat(depot_ind,[skip_nums_extra_scens indices_skips])))'...
    + numV_D*dist_mat(dump_ind,depot_ind));
% Daily operation cost
capital_cost =  skip_add_cost*sum(add_bins_vect) + truck_cost*numV;
total_op_cost = period_op_cost.*numV_D; %sum(ot)
Objective = sum(distance_diff) + sum(total_op_cost) + capital_cost;

%% Set options for YALMIP and solver - Solve the problem
% CUTSDP or gurobi
options =   sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1,'gurobi.MIPGap',0.05);
disp('optimize being called');
sol =       optimize(Constraints,Objective,options);
disp('optimized finished');
toc
%% Analyze error flags & Get the solution
if sol.problem == 0  
    disp('Finished running');
%% Solution processing

    %% Plot results !!! the y axis is not the bin number but the nth smallest scenario index !!!
    close('all')

    figure()
    current_vals = value(xit);
    markers = ['.','*','x','+'];
    legend_group = zeros(1,4);
     
    for i = 1:length(indices_skips)
        extra_period = round(value(yis{i})'*(1./scen_cells{i,3}));
        current_val = round(value(xit(:,i)));
        x = linspace(0.25,T-0.25,T*P);
        y = round(current_val'*i);
       
        current_scatter = scatter(x(y~=0),nonzeros(y),markers(extra_period),'k');%filling_rates(i,1)*50/0.5);
        if ~legend_group(extra_period)
            legend_group(extra_period) = current_scatter;
        end
        hold on
    end
    
    for i = 1:length(skip_nums_extra_scens)
        current_val = round(value(xit_extra(:,i)));
        x = linspace(0.25,T-0.25,T*P);
        y = round(current_val'*(i+length(indices_skips)));
        current_scatter = scatter(x(y~=0),nonzeros(y),markers(extra_period),'k');%filling_rates(i,1)*50/0.5);
        hold on
    end
    
    
    % ADD EXTRA SKIPS HERE
    xlim([0 7])
%     lgd=legend(legend_group(2:4),{'2','3','4'});
%     lgd.Title.String = 'week periodicity';
    xticks([0 1 2 3 4 5 6 7])
    xticklabels({'Mon','Tue','Wed','Thu','Fri','Sat','Sun'})
    ax = gca;
    ax.LineWidth = 0.5;
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ax.GridLineStyle = '-';
    ax.GridColor = 'k';
    ax.GridAlpha = 1;
    ax.YAxis.MinorTickValues = 0:70;
    grid minor

    ylabel('Skip number')
    xlabel('Day')
    grid()
else
    disp('Hmm, something went wrong!');
    sol.info;
    yalmiperror(sol.problem)
end

toc;

%% Analysis
skip_nums_extra_scens.*(sum(value(xit_extra)) ~= 0)
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