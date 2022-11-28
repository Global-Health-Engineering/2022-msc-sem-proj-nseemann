%% Parameters
clear
params.optiParams.period_op_cost = 5400;
params.optiParams.truck_cost = 0;
params.optiParams.km_cost = 384;
params.optiParams.skip_add_cost = 0;
params.optiParams.max_add_bins = 5;
params.optiParams.period_t_max = 4;
params.optiParams.numV = 3;

params.preParams.speed_avg = 30;
params.preParams.dump_ind = 54;
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

clearvars -except params
%%
%load(fullfile(pwd,'test_base','22-11-21_19-44.mat'));
run('full_size_intra_day.m');
%% Save all
ans = output.solver_output.solveroutput.result.objval
mkdir(fullfile('results',[datestr(datetime,'yy-mm-dd') '_output']));
save(fullfile(pwd,'results', [datestr(datetime,'yy-mm-dd') '_output'], [datestr(datetime,'yy-mm-dd_HH-MM') '_output']), 'output', 'params');