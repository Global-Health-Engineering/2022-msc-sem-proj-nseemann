%%%ETH ZURICH 
%%%Course: SP

clear;
tic;
%% Distance matrix
cd 'D:\Semester_project_git\2022-msc-sem-proj-nseemann\src\MATLAB_optim'

%% Import distances

fileID = fopen('../../data/interm_data/dist_matrix_ID_filtered.csv');
formatSpec = '%s';
N = 3;
C_text = textscan(fileID,formatSpec,N,'Delimiter',',');

C = textscan(fileID,'%d %d %f','Delimiter',',');

fclose(fileID);



%% Distances to matrix
dist_mat = zeros(length(max(max(C{1},max(C{2})))));
for i = 1:length(C{1})
    dist_mat(C{1}(i),C{2}(i)) = C{3}(i);
end

%% Parameters
numBins = 3;

filling_rates = [0.09, 0.5, 0.3];%rand(1,numBins)/2;

%Reduced weekly scenarios, only for two visits per week and period 1 week

[scens,scensgaps]=create_scens();

scensize = size(scens);

scenlen = scensize(1);
T = scensize(2);

dump_ind = 54; % Mzedi dump
compost_ind = 55; % Limbe composting
depot_ind = 56; % Storage city center
%    d  du s1 s2 s3
dist =[0 4 1 2 3;
    4 0 1 3 4;
    1 1 0 2 9;
    2 3 2 0 5;
    3 4 9 5 0];
day_op_cost = 20;
%% Decision variables
xit = binvar(T,numBins,'full'); %if is operating on day t
yis = binvar(scenlen, numBins,'full'); %if is assigned scenario s
fit = binvar(T, numBins,'full'); %if first on day t (can have multiple)
ot = binvar(T,1); % if operating on day t
%% Constraints


% Equality constraints
assign_1 = sum(yis) == 1; % Assigns one scenario to each skip
day_flow = xit'-yis'*scens == 0; %Bins emptied once on days of schedule

% Inequality constraints
day_op = ot >= sum(xit,2)/numBins; %Counts operation days
first = sum(fit,2) == ot; %Makes it so at least one is first

fill_min = scensgaps*yis >= filling_rates; %All bins at least on schedule

% Overall constraints
Constraints = [assign_1 day_flow fill_min day_op first];
%% Objective function

distance_diff = fit*((-1*dist(dump_ind,3:dump_ind+numBins) + dist(depot_ind,3:dump_ind+numBins)))';
Objective = sum(sum(xit)) + sum(ot)*day_op_cost + sum(distance_diff);


%% Set options for YALMIP and solver - Solve the problem
options =   sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1);
sol =       optimize(Constraints,Objective,options);

fprintf('Total system cost: = %d', sol.solveroutput.result.objval);

%% Analyze error flags & GET the solution
if sol.problem ==0  
    disp('All good');
else
    disp('Hmm, something went wrong!');
    sol.info;
    yalmiperror(sol.problem)
end

%%

%% Functions
% Scenarios

function [scens,sensgaps] = create_scens()
    T=28;
    sensgaps = [1/2];
     scens = repmat([1 1 1 1 1 1 0],1,4);
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
