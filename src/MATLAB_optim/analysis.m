clear
%% Which extra skips:
load(fullfile(pwd,'results','22-11-28_output','22-11-28_18-11_output.mat'));

add_bins_vect_1_11 = output.results.add_bins_vect;

load(fullfile(pwd,'results','22-11-29_output','22-11-29_13-54_output.mat'));

add_bins_vect_2_11 = output.results.add_bins_vect;

[add_bins_vect_1_11' add_bins_vect_2_11' output.process.filling_rates(:,1)]

%% Number of rounds depending on added skips
%%Costs per day depending on added skips
%%Costs per round depending on added

total_rounds = sum(output.optiVars.numV_D)
dist_cost_per_round = sum(output.results.distance_cost)./total_rounds
%% Number of possible extra skips from scenarios
possible_extra = 0;
for i = 1:length(output.process.scen_cells)
    possible_extra = possible_extra + max(output.process.scen_cells{i,2});
end
possible_extra
%% Finding added skips
% nums of skips added
skips_added_indices = output.process.skip_nums_extra_scens(find(sum(output.optiVars.xit_extra)))

% filling rate of respective skips
output.process.filling_rates(output.process.skip_nums_extra_scens(find(sum(output.optiVars.xit_extra))))

% distance dump-skip
(output.process.dist_mat(params.preParams.dump_ind,skips_added_indices)' + output.process.dist_mat(skips_added_indices,params.preParams.dump_ind))/2
%% Comparing scen_cell collection intensity depending on scenario and added skips
i=16
[sum(output.process.scens(output.process.scen_cells{i,1},:),2) output.process.scen_cells{i,2} output.process.filling_rates(i,1)./(1+output.process.scen_cells{i,2})]

%%
for i = 1:length(output.process.scen_cells)
    current_array = output.process.scen_cells{i,2};
    for j = 1:length(current_array)
        if current_array(j) == 2
            disp('is2');
        end
    end
end