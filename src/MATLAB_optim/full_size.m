%%%ETH ZURICH 
%%%Course: SP

clear;
tic;

cd 'D:\Semester_project_git\2022-msc-sem-proj-nseemann\src\MATLAB_optim'

%% Variables

% Cost of operating for 1 day
day_op_cost = 20;


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

%% Filtering distances and isolating dump, composting and storage facility



indices_util = sort([dump_ind,compost_ind,depot_ind]);
indices_skips = [1:indices_util(1)-1,indices_util(1)+1:indices_util(2)-1,indices_util(2)+1:indices_util(3)-1,indices_util(3)+1:max(max(C{1},max(C{2})))];

dist_mat_no_util = dist_mat([1:indices_util(1)-1,indices_util(1)+1:indices_util(2)-1,indices_util(2)+1:indices_util(3)-1,indices_util(3)+1:max(max(C{1},max(C{2})))],[1:indices_util(1)-1,indices_util(1)+1:indices_util(2)-1,indices_util(2)+1:indices_util(3)-1,indices_util(3)+1:max(max(C{1},max(C{2})))]);
%% Parameters

numBins = length(indices_skips);

% Random but constant rates
%load('filling_rates_13_10_22.mat')
filling_rates = rand(1,numBins)/2;
%% Load scenarios and relevant variables

[scens,scensgaps]=create_scens();

scensize = size(scens);

scenlen = scensize(1);
T = scensize(2);

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
first = sum(fit, 2) == ot; %Makes it so at least one is first
fit_unity = fit <= xit;

fill_min = scensgaps*yis >= filling_rates; %All bins at least on schedule

% Overall constraints
Constraints = [assign_1 day_flow fill_min day_op first fit_unity];
%% Objective function

distance_diff = xit*(dist_mat(dump_ind,indices_skips))' + ((dist_mat(dump_ind,indices_skips))*xit')' + fit*((-1*dist_mat(dump_ind,indices_skips) + dist_mat(depot_ind,indices_skips)))';
Objective = sum(distance_diff) + 20*sum(ot);%sum(sum(xit)) + sum(ot)*day_op_cost + sum(distance_diff);

%% Set options for YALMIP and solver - Solve the problem
options =   sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1);
sol =       optimize(Constraints,Objective,options);

scen_assignment = zeros(2, length(indices_util) + length(indices_skips));
b = 1;
value_scen_assignment = value(yis);
for i = 1:length(scen_assignment)
    if ismember(i,indices_util)
        scen_assignment(1,i) = NaN;
        scen_assignment(2,i) = NaN;
    elseif ismember(i, indices_skips)
        scen_assignment(1,i) = find(value_scen_assignment(:,b));
        scen_assignment(2,i) = filling_rates(b);
        b = b+1;
    end
end

disp(scen_assignment)

%fprintf('Total system cost: = %d', sol.solveroutput.result.objval);

%% Analyze error flags & GET the solution
if sol.problem == 0  
    disp('Finished running');
else
    disp('Hmm, something went wrong!');
    sol.info;
    yalmiperror(sol.problem)
end

%% Plot results
current_vals = value(xit);
current_vals_scens = [current_vals' rmmissing(scen_assignment(1,:))'];
for i = 1:length(indices_skips)
    %bin_ind = indices_skips(i);
    %current_val = current_vals_scens(bin_ind,:);
    [~,ind] = min(current_vals_scens(:,29));
    %disp(ind)
    current_val = current_vals_scens(ind,1:28);
    x = linspace(1,28,28);
    y = current_val'*i;
    scatter(x(y~=0),nonzeros(y));
    current_vals_scens(ind,:) = [];
    %disp(current_vals_scens(:,29))
    hold on
end

toc
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
