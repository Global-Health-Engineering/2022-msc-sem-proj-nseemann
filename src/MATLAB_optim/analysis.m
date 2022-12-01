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

%%
possible_extra = 0;
for i = 1:length(output.process.scen_cells)
    possible_extra = possible_extra + max(output.process.scen_cells{i,2});
end
possible_extra
%%
indices_05 = find(output.process.filling_rates(:,1) == 0.5);
indices_05_1 = indices_05(7);
indices_0_skips = output.process.scen_cells{indices_05_1,1}(output.process.scen_cells{indices_05_1,2} == 0);
size(indices_0_skips)
output.process.scens(indices_0_skips,:)

indices_1_skips = output.process.scen_cells{indices_05_1,1}(output.process.scen_cells{indices_05_1,2} == 1);
size(indices_1_skips)
output.process.scens(indices_1_skips,:)
output.process.scens(output.optiVars.yis{indices_05_1}'*output.process.scen_cells{indices_05_1,1},:)
%%
for i = 1:length(output.process.scen_cells)
    current_array = output.process.scen_cells{i,2};
    for j = 1:length(current_array)
        if current_array(j) == 2
            disp('is2');
        end
    end
end