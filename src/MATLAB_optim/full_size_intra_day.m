%%%ETH ZURICH 
%%%Course: SP
close all
clear; 
tic;

% Change to current running directory:
%cd 'D:\Semester_project_git\2022-msc-sem-proj-nseemann\src\MATLAB_optim'

%% Variables

% Cost of operating one truck for 1 day (assume labour costs)
period_op_cost = 4.5;

% Average speed 30km/h
speed_avg = 30; 

% Capital cost of truck / maximum number of trucks
truck_cost = 0;
max_truck = 2;

% Distance cost
km_cost = 0.5;

% Capital cost of buying one skip
skip_add_cost = 0;% 300;
max_add_bins = 0;

% Indices other than skips
dump_ind = 54; % Mzedi dump
compost_ind = 55; % Limbe composting
depot_ind = 56; % Storage city center


%% Import distances

%Distance matrix
fileID = fopen('../../data/interm_data/dist_matrix_ID_filtered.csv');
textscan(fileID,'%s',3,'Delimiter',',');
C = textscan(fileID,'%d %d %f','Delimiter',',');
fclose(fileID);

%Convert into a MATLAB matrix
dist_mat = zeros(length(max(max(C{1},max(C{2})))));
for i = 1:length(C{1})
    dist_mat(C{1}(i),C{2}(i)) = C{3}(i);
end

%Time matrix
t_mat = dist_mat/speed_avg;

%% Filtering distances and isolating dump, composting and storage facility
indices_util = sort([dump_ind,compost_ind,depot_ind]);
indices_skips = [1:indices_util(1)-1,indices_util(1)+1:indices_util(2)-1,indices_util(2)+1:indices_util(3)-1,indices_util(3)+1:max(max(C{1},max(C{2})))];

%dist_mat_no_util = dist_mat([1:indices_util(1)-1,indices_util(1)+1:indices_util(2)-1,indices_util(2)+1:indices_util(3)-1,indices_util(3)+1:max(max(C{1},max(C{2})))],[1:indices_util(1)-1,indices_util(1)+1:indices_util(2)-1,indices_util(2)+1:indices_util(3)-1,indices_util(3)+1:max(max(C{1},max(C{2})))]);
%% Parameters

numBins = length(indices_skips);

% Random but constant rates
load('filling_rates_13_10_22.mat')
%figure()
% Create randoms:
%filling_rates = rand(1,numBins)/2;
%histogram(filling_rates)
%ylabel('number of skips')
%xlabel('filling rate [/day]')
%% Load scenarios and relevant variables

% Load scenarios and largest gap in scenario
[scens,mult]=create_scens();
divis = 1./mult;
% Get period and number of scenarios
scensize = size(scens);
% Time days
T = 7;
% Number of collections periods per day
P = 2;

% maximum number of additional bins at each bin
set_add_bins = 2;

% additional bins only considered for more or equal than per_week_consider services
% under current service
per_week_consider = 3;

% Makes uniform filling rates over entire period
load('filling_rates_13_10_22.mat');
filling_rates = repmat(filling_rates',1,T*P);
filling_rates(1,:) = [repmat([0.66],1, 14)];% repmat([0.3], 1,3)];
scennum = scensize(1);

%% Narrow relevant scenarios
scen_cells = {};
%if a scenarios is suitable, add the scenario vector to skip scenario cell
excessive = [];
excess_bins = []; 
for i = 1:numBins
    scen_cells{i,1} = [];
    scen_cells{i,2} = [];
    scen_cells{i,3} = [];
    current_rates_doubled = [filling_rates(i,:) filling_rates(i,:)]; %doubled for extra weekly
    for s = 1:size(scens,1)
        current_scen_doubled = [scens(s,:) scens(s,:)]; %doubled for extra weekly
        previous = 1;
        exceed = 0;
        sum_period = [];
        k = 1;
        for d = find(current_scen_doubled)
            sum_period = [sum_period sum(current_rates_doubled(previous:d)*mult(s))/P];
            if sum_period(end) >= 1
                exceed = 1;
                break
            end
            k = k+1;
            previous = d + 1;
        end 
        if exceed ~= 1
            scen_cells{i,1} = [scen_cells{i,1}; s]; %scens(s,:) <-- old entire scen
            scen_cells{i,2} = [scen_cells{i,2}; 0];
            scen_cells{i,3} = [scen_cells{i,3}; divis(s)];
        end
    end
    for add_bins = 1:set_add_bins
        %Scens with more than per_week_consider in a week
        if min(sum(scens(scen_cells{i,1},:),2).*scen_cells{i,3}) >= per_week_consider 
            current_rates_doubled = [filling_rates(i,:)./(1+add_bins) filling_rates(i,:)./(1+add_bins)]; %doubled for extra weekly
            excess_bins = [excess_bins; i add_bins filling_rates(i,1)./(1+add_bins)];
            for s = 1:size(scens,1)
                current_scen_doubled = [scens(s,:) scens(s,:)]; %doubled for extra weekly
                previous = 1;
                exceed = 0;
                sum_period = [];
                k = 1;
                for d = find(current_scen_doubled)
                    sum_period = [sum_period sum(current_rates_doubled(previous:d)*mult(s))/P];
                    if sum_period(end) >= 1
                        exceed = 1;
                        break
                    end
                    k = k+1;
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

%% Decision variables
xit = binvar(T*P,numBins,'full'); %if is operating on day t at period p

yis = {};
for l = 1:size(scen_cells,1)
    yis{l} = binvar(size(scen_cells{l,1},1), 1,'full'); %if l is assigned scenario in subset
end

fit = binvar(T*P, numBins,'full'); %if first on day t at p
ot = binvar(T*P,1,'full'); % if operating on day t at p
numV_D = intvar(T*P,1,'full'); % Number of vehicles in operation on day

numV = intvar(1); % Total number of vehicles

%% Constraints
% Equality constraints

assign_1 = []; %Only one scenario within allowed scenario per skip
for l = 1:size(scen_cells,1)
    assign_1 = [assign_1 sum(yis{l})==1];
end

day_flow = []; %Skips emptied only on periods assigned by scenario
for l = 1:size(scen_cells,1)
    day_flow = [day_flow xit(:,l)'-yis{l}'*scens(scen_cells{l,1},:)==0];
end

period_t_max = 3;


% Inequality constraints

% Ensure period time is not exceeded
period_time = xit*(t_mat(dump_ind,indices_skips))' + ((t_mat(dump_ind,indices_skips))*xit')' + fit*((-1*t_mat(dump_ind,indices_skips) + t_mat(depot_ind,indices_skips)))' <= repmat(period_t_max,T*P,1).*numV_D;

period_op = ot >= sum(xit,2)/numBins; %Counts operation days (periods)
first = sum(fit, 2) == numV_D; %Makes it so numV_D skip is first in loop on each period
fit_unity = fit <= xit; %Prevents a first in loop when a skip is not operating on that period

%fill_min = scensgaps*yis >= filling_rates; %All bins at least on schedule

numV_const = 1<= numV <= max_truck;
numV_D_numV = zeros(T*P,1) <= numV_D <= repmat(numV,T*P,1);

total_added_skips = 0;
for l = 1:size(scen_cells,1)
   total_added_skips = total_added_skips + yis{l}'*scen_cells{l,2};
end

add_skip_max = total_added_skips <= max_add_bins;
% Overall constraints
Constraints = [assign_1 numV_D_numV day_flow period_op first fit_unity period_time numV_const add_skip_max];
%fill_min
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

% Distance travelled on the rounds
distance_diff = 20*sum(divis_vect) +  xit*(dist_mat(dump_ind,indices_skips))' + ((dist_mat(dump_ind,indices_skips))*xit')' + fit*((-1*dist_mat(dump_ind,indices_skips) + dist_mat(depot_ind,indices_skips)))';
% Daily operation cost
capital_cost =  skip_add_cost*sum(add_bins_vect) + truck_cost*numV;
total_op_cost = period_op_cost*sum(ot);
Objective = sum(distance_diff) + total_op_cost + capital_cost;

%% Set options for YALMIP and solver - Solve the problem
% CUTSDP or gurobi
options =   sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1,'gurobi.MIPGap',0.005);
sol =       optimize(Constraints,Objective,options);
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
        extra_period = value(yis{i}).'*(1./scen_cells{i,3});
        current_val = value(xit(:,i));
        x = linspace(0.25,T-0.25,T*P);
        y = current_val'*i;
       
        current_scatter = scatter(x(y~=0),nonzeros(y),markers(extra_period),'k');%filling_rates(i,1)*50/0.5);
        if ~legend_group(extra_period)
            legend_group(extra_period) = current_scatter;
        end
        hold on
    end
    xlim([0 7])
    legend(legend_group,{'1','2','3','4'});
    xticks([0 1 2 3 4 5 6 7])
    xticklabels({'Mon','Tue','Wed','Thu','Fri','Sat','Sun'})
    ax = gca;
    ax.LineWidth = 0.5;
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ax.GridLineStyle = '-';
    ax.GridColor = 'k';
    ax.GridAlpha = 1;
    ax.YAxis.MinorTickValues = 0:60;
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
%%
% [scens,mult] = create_scens();
% size(scens)
%% Scenarios 1x14 double period one weekly with week multiplier/cost divisor
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