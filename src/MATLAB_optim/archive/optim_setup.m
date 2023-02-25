%%%ETH ZURICH 
%%%Course: SP

clear;
tic;

%% Parameters
numBins = 3;
T = 7;

filling_rates = [0.09, 0.5, 0.3];%rand(1,numBins)/2;

%Reduced weekly scenarios, only for two visits per week and period 1 week
scenarios= [[1 0 1 0 1 0 0];[1 1 1 1 1 1 1];[0 1 0 0 0 0 0];[1 0 0 0 0 0 0]];

[scens,scensgaps]=create_scens();

scensize = size(scenarios);

scenlen = scensize(1);

gs = [1/3 1 1/7 1/7];
%    d  du s1 s2 s3
d =[0 4 1 2 3;
    4 0 1 3 4;
    1 1 0 2 9;
    2 3 2 0 5;
    3 4 9 5 0];

%% Decision variables
xit = binvar(T,numBins,'full'); 
yis = binvar(scenlen, numBins,'full');
fit = binvar(T, numBins,'full');

%% Constraints


% Equality constraints
assign_1 = sum(yis) == 1;
day_flow = xit'-yis'*scenarios == 0;
% Inequality quality

fill_min = gs*yis >= filling_rates;


% Overall constraints
Constraints = [assign_1 day_flow fill_min];
%% Objective function

Objective = sum(sum(xit));


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

%% fit*((-1*dist_mat(dump_ind,indices_skips) 

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
         sensgaps=[sensgaps 1/2];
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
