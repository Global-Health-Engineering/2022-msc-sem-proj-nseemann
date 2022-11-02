%%%ETH ZURICH 
%%%Course: SP
close all
clear; 
tic;

% Change to current running directory:
%cd 'D:\Semester_project_git\2022-msc-sem-proj-nseemann\src\MATLAB_optim'

%% Variables

% Cost of operating one truck for 1 day
day_op_cost = 20;

% Average speed 30km/h
speed_avg = 30; 

% Capital cost of truck
truck_cost = 1000;

% Operation cost day per truck
truck_cost_D = 50;

% Capital cost of buying one skip
skip_cost = 200;

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
[scens,scensgaps]=create_scens();

% Get period and number of scenarios
scensize = size(scens);
% Time days
T = 28;
% Number of collections per day
P = 2;
scennum = scensize(1);

%% Decision variables
xit = binvar(T*P,numBins,'full'); %if is operating on day t at period p
yis = binvar(scennum, numBins,'full'); %if is assigned scenario s
fit = binvar(T*P, numBins,'full'); %if first on day t at p
ot = binvar(T*P,1,'full'); % if operating on day t at p
numV_D = intvar(T*P,1,'full'); % Number of vehicles in operation on day

numV = intvar(1); % Total number of vehicles

% Constraints
% Equality constraints
assign_1 = sum(yis) == 1; % Assigns one scenario to each skip
day_flow = xit'-yis'*scens == 0; %Bins emptied once on period of schedule


day_t_max = 8;

day_time = xit*(t_mat(dump_ind,indices_skips))' + ((t_mat(dump_ind,indices_skips))*xit')' + fit*((-1*t_mat(dump_ind,indices_skips) + t_mat(depot_ind,indices_skips)))' <= repmat(day_t_max,T*P,1).*numV_D;
% Inequality constraints
day_op = ot >= sum(xit,2)/numBins; %Counts operation days (periods)
first = sum(fit, 2) == ot.*numV_D; %Makes it so numV_D skip is first in loop on each day
fit_unity = fit <= xit; %Prevents a first in loop when a skip is not operating on that day

%fill_min = scensgaps*yis >= filling_rates; %All bins at least on schedule

numV_min = numV >= 1;
numV_D_numV = zeros(T*P,1) <= numV_D <= repmat(numV,T*P,1);

% Overall constraints
Constraints = [assign_1 numV_D_numV day_flow day_op first fit_unity day_time numV_min ];
%fill_min
%% Objective function
% Distance travelled on the rounds
distance_diff = xit*(dist_mat(dump_ind,indices_skips))' + ((dist_mat(dump_ind,indices_skips))*xit')' + fit*((-1*dist_mat(dump_ind,indices_skips) + dist_mat(depot_ind,indices_skips)))';
% Daily operation cost
capital_cost =  truck_cost*numV;
total_op_cost = day_op_cost*sum(ot) + sum(truck_cost_D*numV_D);
Objective = sum(distance_diff) + total_op_cost + capital_cost;

%% Set options for YALMIP and solver - Solve the problem
options =   sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1);
sol =       optimize(Constraints,Objective,options);


%% Analyze error flags & Get the solution
if sol.problem == 0  
    disp('Finished running');
        %% Solution processing
    % Scen assignment to output matrix coupling the filling rate with the
    % scenario number
    scen_assignment = zeros(2, length(indices_util) + length(indices_skips));
    b = 1;
    value_scen_assignment = value(yis);
    for i = 1:length(scen_assignment)
        if ismember(i,indices_util)
            scen_assignment(1,i) = NaN;
            scen_assignment(2,i) = NaN;
        elseif ismember(i, indices_skips)
            scen_assignment(1,i) = find(value_scen_assignment(:,b),1);
            scen_assignment(2,i) = filling_rates(b);
            b = b+1;
        end
    end

    %% Plot results !!! the y axis is not the bin number but the nth smallest scenario index !!!
    figure()
    current_vals = value(xit);
    current_vals_scens = [current_vals' rmmissing(scen_assignment(1,:))'];
    disp(current_vals_scens)
    for i = 1:length(indices_skips)
        %bin_ind = indices_skips(i);
        %current_val = current_vals_scens(bin_ind,:);
        [~,ind] = min(current_vals_scens(:,57)); 
        %disp(ind)
        current_val = current_vals_scens(ind,1:T*P);
        %disp([ 'current val: ' num2str(scen_assignment(1,ind))])
        x = linspace(1,T*P,T*P);
        y = current_val'*i;
        scatter(x(y~=0),nonzeros(y));
        current_vals_scens(ind,:) = [];
        %disp(current_vals_scens(:,29))
        hold on
    end


    ylabel('Skip number')
    xlabel('Day')
    grid()
else
    disp('Hmm, something went wrong!');
    sol.info;
    yalmiperror(sol.problem)
end



toc
%% Scenarios 1x28 basics

function [scens,sensgaps] = create_scens()
      scens = repmat([1 1 1 1 1 1 1 1 1 1 1 1 0 0],1,4);
      scens = [scens;repmat([1 0 0 0 0 0 0 1 0 0 0 0 0 0],1,4)];
      scens = [scens;repmat([ 0 1 0 0 0 0 0 0 1 0 0 0 0 0],1,4)];
      scens = [scens;repmat([ 0 0 1 0 0 0 0 0 0 1 0 0 0 0],1,4)];
      scens = [scens;repmat([ 0 0 0 1 0 0 0 0 0 0 1 0 0 0],1,4)];
      scens = [scens;repmat([ 0 0 0 0 1 0 0 0 0 0 0 1 0 0],1,4)];
      sensgaps = 4;
%     T=28;
%     % Everyday
%     sensgaps = [1/2];
%     scens = repmat([1 1 1 1 1 1 0],1,4);
%     
%      %once a month:
%      N_f = 28;
%      for i = 1:N_f
%          if rem(i,7) ~= 0
%              sensgaps = [sensgaps 1/N_f];
%             scens = [scens; [ repmat([0],1,i-1) 1 repmat([0], 1,N_f-i)]];
%          end
%      end
%      %2/month
%      N_f = 14;
%      for i = 1:N_f
%          if rem(i,7) ~= 0
%              sensgaps = [sensgaps 1/N_f];
%             scens = [scens; repmat([repmat([0],1,i-1) 1 repmat([0],1,N_f-i)],1,T/N_f)];
%          end
%      end 
%      
%      %4/month
%      N_f = 7;
%      for i = 1:N_f
%          if rem(i,7) ~= 0
%              sensgaps = [sensgaps 1/N_f];
%             scens = [scens; repmat([repmat([0],1,i-1) 1 repmat([0],1,N_f-i)],1,T/N_f)];
%          end
%      end 
%      
%      %2/week
%      basis = [1 0 0 1 0 0 0;0 1 0 0 0 1 0; 1 0 0 0 1 0 0;  0 0 1 0 0 1 0; 0 1 0 0 1 0 0 ];
%      scens = [scens; repmat(basis,1,4)];
%      for i = 1:5
%          sensgaps=[sensgaps 1/3];
%      end
%      
%      %3/week
%      basis = [1 0 1 0 1 0 0; 0 1 0 1 0 1 0; 1 0 0 1 0 1 0; 1 0 1 0 0 1 0];
%      scens = [scens; repmat(basis,1,4)];
%      for i = 1:4
%          sensgaps=[sensgaps 1/3];
%      end
%      
%      % 4/week
%      basis = [1 0 1 0 1 1 0; 1 1 0 1 0 1 0; 1 0 1 1 0 1 0];
%      scens = [scens; repmat(basis,1,4)];
%      for i = 1:3
%          sensgaps=[sensgaps 1/2];
%      end
%      
%      % 5/week
%      basis = [1 0 1 1 1 1 0; 1 1 0 1 1 1 0; 1 1 1 0 1 1 0; 1 1 1 1 0 1 0];
%      scens = [scens; repmat(basis,1,4)];
%      for i = 1:4
%          sensgaps=[sensgaps 1/2];
%      end
%      % 6/week : already at the top!
     
end 

%% Scenarios 2x28 limited range

function [scens,sensgaps] = create_scens_enhanced()
    %protocol:
    %create cell array with 2x28 entries for each scenarios
    %create array with N_scen scenarios chatacteristics
    %take into account Sunday might have different filling rate than during
    %the week
    T=28;
    
    % Everyday
    sensgaps = [1/2];
    scens = repmat([1 1 1 1 1 1 0; 0 0 0 0 0 0 0],1,4);
    scens = [scensrepmat([1 1 1 1 1 1 0; 0 0 0 0 0 0 0],1,4)];
    
     %once a month:
     N_f = 28;
     
     for i = 1:N_f
         if rem(i,7) ~= 0
             sensgaps = [sensgaps 1/N_f];
            scens = [scens; [ repmat([0],1,i-1) 1 repmat([0], 1,N_f-i)]];
         end
     end
     %2/month
     N_f = 14;
     for i = 1:N_f
         if rem(i,7) ~= 0
             sensgaps = [sensgaps 1/N_f];
            scens = [scens; repmat([repmat([0],1,i-1) 1 repmat([0],1,N_f-i)],1,T/N_f)];
         end
     end 
     
     %4/month
     N_f = 7;
     for i = 1:N_f
         if rem(i,7) ~= 0
             sensgaps = [sensgaps 1/N_f];
            scens = [scens; repmat([repmat([0],1,i-1) 1 repmat([0],1,N_f-i)],1,T/N_f)];
         end
     end 
     
     %2/week
     basis = [1 0 0 1 0 0 0;0 1 0 0 0 1 0; 1 0 0 0 1 0 0;  0 0 1 0 0 1 0; 0 1 0 0 1 0 0 ];
     scens = [scens; repmat(basis,1,4)];
     for i = 1:5
         sensgaps=[sensgaps 1/3];
     end
     
     %3/week
     basis = [1 0 1 0 1 0 0; 0 1 0 1 0 1 0; 1 0 0 1 0 1 0; 1 0 1 0 0 1 0];
     scens = [scens; repmat(basis,1,4)];
     for i = 1:4
         sensgaps=[sensgaps 1/3];
     end
     
     % 4/week
     basis = [1 0 1 0 1 1 0; 1 1 0 1 0 1 0; 1 0 1 1 0 1 0];
     scens = [scens; repmat(basis,1,4)];
     for i = 1:3
         sensgaps=[sensgaps 1/2];
     end
     
     % 5/week
     basis = [1 0 1 1 1 1 0; 1 1 0 1 1 1 0; 1 1 1 0 1 1 0; 1 1 1 1 0 1 0];
     scens = [scens; repmat(basis,1,4)];
     for i = 1:4
         sensgaps=[sensgaps 1/2];
     end
     % 6/week : already at the top!
     
end 