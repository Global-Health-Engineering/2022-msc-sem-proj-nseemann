clear
%% Which extra skips:
load(fullfile(pwd,'results','22-11-28_output','22-11-28_18-11_output.mat'));

add_bins_vect_1_11 = output.results.add_bins_vect;

load(fullfile(pwd,'results','22-11-29_output','22-11-29_13-54_output.mat'));

add_bins_vect_2_11 = output.results.add_bins_vect;

[add_bins_vect_1_11' add_bins_vect_2_11' output.process.filling_rates(:,1)]

%%
load(fullfile(pwd,'results','22-11-28_output','22-11-28_18-11_output.mat'));
indices_06 = find(output.process.filling_rates(:,1) > 0.6);
indices_06_1 = indices_06(1);
indices_0_skips = output.process.scen_cells{indices_06_1,1}(output.process.scen_cells{indices_06_1,2} == 0);
size(indices_0_skips)
output.process.scens(indices_0_skips,:)

indices_1_skips = output.process.scen_cells{indices_06_1,1}(output.process.scen_cells{indices_06_1,2} == 1);
size(indices_1_skips)
output.process.scens(indices_1_skips,:)

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
%%
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