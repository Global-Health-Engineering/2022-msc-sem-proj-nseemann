%% Parameters
clear
params.optiParams.period_op_cost = 5400;
params.optiParams.truck_cost = 0;
params.optiParams.km_cost = 384;
params.optiParams.skip_add_cost = 0;
params.optiParams.max_add_bins = 0;
params.optiParams.period_t_max = 4;
params.optiParams.numV = 2;

params.preParams.speed_avg = 30;
params.preParams.dump_ind = 54; %normally 54, but can replace with compost = 55
params.preParams.compost_ind = 55;
params.preParams.depot_ind = 56;
params.preParams.T = 7;
params.preParams.P = 2;
params.preParams.underfull_threshold = 0.6;
params.preParams.set_add_bins = 2;
params.preParams.per_week_consider = 3;
params.preParams.unified_salary = 1;
params.preParams.nanreplace = 0.2;

params.optiOptions = sdpsettings('verbose',1,'solver','gurobi','savesolveroutput',1,'gurobi.MIPGap',0.01);

params_basis = params;
clearvars -except params params_basis 

%% 
clearvars -except params_basis
experiment = zeros(1,7); %numV,set_add,max_add,possible_add,added,objective,gap
params = params_basis;

params.optiParams.numV = 4;
params.optiParams.period_t_max = 4;
params.preParams.set_add_bins = 1;
vary_max_add_bins = 2;

%table_parameters = ["numV";"set_add_bins";"possible_extra";"feas","added";"objective";"gap"];

ID = cell(length(vary_max_add_bins),1);
numV = zeros(length(vary_max_add_bins),1);
set_add_bins = zeros(length(vary_max_add_bins),1);
max_add = zeros(length(vary_max_add_bins),1);
possible_extra = zeros(length(vary_max_add_bins),1);
feas = zeros(length(vary_max_add_bins),1);
added = zeros(length(vary_max_add_bins),1);
objective = zeros(length(vary_max_add_bins),1);
gap = zeros(length(vary_max_add_bins),1);

b = 1;
for i = vary_max_add_bins
    params.optiParams.max_add_bins = i;
    run('full_size_intra_day.m');
    current_datetime = [datestr(datetime,'yy-mm-dd')];
    current_datetime_ext = [datestr(datetime,'yy-mm-dd_HH-MM') num2str(b)];
    ID{b} = [current_datetime_ext '_output'];
    numV(b) = params.optiParams.numV;
    set_add_bins(b) = params.preParams.set_add_bins;
    max_add(b) = params.optiParams.max_add_bins;
    possible_extra(b) = output.process.possible_extra;
    feas(b) = output.feasible;
    if output.feasible == 1
        added(b) = sum(output.results.add_bins_vect);
        objective(b) = output.results.Objective;
        gap(b) = output.bound_gap;
    end
    b = b+1;
    mkdir(fullfile('results',[current_datetime '_output']));
    save(fullfile(pwd,'results', [current_datetime '_output'], [current_datetime_ext '_output']), 'output', 'params'); 
end

table_exp = table(numV,set_add_bins,max_add,feas,added,objective,gap,'RowName',ID);
writetable(table_exp,'2_add_per_week_3');
%%
%load(fullfile(pwd,'test_base','22-11-21_19-44.mat'));
%run('full_size_intra_day.m');
%% Save all
%ans = output.solver_output.solveroutput.result.objval
% mkdir(fullfile('results',[datestr(datetime,'yy-mm-dd') '_output']));
% save(fullfile(pwd,'results', [datestr(datetime,'yy-mm-dd') '_output'], [datestr(datetime,'yy-mm-dd_HH-MM') '_output']), 'output', 'params');