%clear
%load(fullfile(pwd,'results','22-11-25_output','22-11-25_22-22_output.mat'));
%% Plot results
%close all
T = params.preParams.T;
P = params.preParams.P;
%%
figure()
t = tiledlayout(1,1);
current_vals = output.optiVars.xit;
markers = ['.','*','x','+'];
legend_group = zeros(1,4);

for i = 1:length(output.process.indices_skips)
    extra_period = round(value(output.optiVars.yis{i})'*(1./output.process.scen_cells{i,3}));
    current_val = round(output.optiVars.xit(:,i));
    x = linspace(0.25,T-0.25,T*P);
    y = round(current_val'*i);

    current_scatter = scatter(x(y~=0),nonzeros(y),markers(extra_period),'k');%filling_rates(i,1)*50/0.5);
    if ~legend_group(extra_period)
        legend_group(extra_period) = current_scatter;
    end
    hold on
end
b=1;
for i = 1:length(output.process.skip_nums_extra_scens)
    current_val = round(output.optiVars.xit_extra(:,i));
    if sum(current_val) ~= 0
        x = linspace(0.25,T-0.25,T*P);
        y = round(current_val'*(b+1+length(output.process.indices_skips)));
        current_scatter = scatter(x(y~=0),nonzeros(y),markers(extra_period),'k');
        b=b+1;
    end
    hold on
end


% ADD EXTRA SKIPS HERE
xlim([0 7])
xticks(sort([linspace(0.25,6.25,7) linspace(0.75,6.75,7)]))
h = arrayfun(@(a,b)xline(a, "-",'LabelHorizontalAlignment','left','LabelVerticalAlignment','top'),[1 2 3 4 5 6 7]);
xticklabels(repmat(["am","pm"],1,7))
yline(length(output.process.indices_skips)+1)
ylabel('Skip')
xlabel('Period of day')
grid()
ax2 = axes(t);
ax2.XAxisLocation = 'top';
ax2.YAxis.Visible  = 'off';
xticks(ax2, 0.5+linspace(0,6,7))
xticklabels(ax2, {'Mon','Tue','Wed','Thu','Fri','Sat','Sun'});
set(ax2.XAxis,'TickDir', 'both');
%ax2.Color = 'none';
ax2.Box = 'off';
set(ax2, 'XMinorGrid', 'off')

xlim(ax2,[0 7])
lgd=legend(legend_group(2:4),{'2','3','4'},'Location', 'northeast');

lgd.Title.String = 'Period';


 %% Daily operation breakdown analysis
%Day breakdown
no_days = [];

for d = 1:T*P
    joined_xit = [output.optiVars.xit output.optiVars.xit_extra];
    no_days = [no_days sum(joined_xit(d,:))];
    %find(joined_xit(d,:))
end
figure()
t = tiledlayout(1,1);
bar(x, [no_days;output.optiVars.numV_D']');
ylim([0,20])
ylabel('Number of skips serviced | Number of vehicles operated')
xlim([0,7])
xticks(sort([linspace(0.25,6.25,7) linspace(0.75,6.75,7)]))

xticklabels(repmat(["am","pm"],1,7))
%xticklabels({'Mon','Tue','Wed','Thu','Fri','Sat','Sun'})
h = arrayfun(@(a,b)xline(a, "-",'LabelHorizontalAlignment','left','LabelVerticalAlignment','top'),[1 2 3 4 5 6 7]);

ax = gca;
ax.LineWidth = 0.5;
grid minor
set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'on')
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = 1;
ax.YAxis.MinorTickValues = 0:18:18;
ax.XAxis.MinorTickValues = [];
set(gca,'MinorGridLineStyle','-')
set(ax,'TickLength', [0.01 0.01]);
legend('Number of skips serviced','Number of vehicles operated')

ax2 = axes(t);

ax2.XAxisLocation = 'top';
ax2.YAxis.Visible  = 'off';
ax2.XAxis.Visible  = 'on';
xticks(ax2, 0.5+linspace(0,6,7))
xticklabels(ax2, {'Mon','Tue','Wed','Thu','Fri','Sat','Sun'});
ax2.Box = 'on';
set(ax2, 'XMinorGrid', 'off', 'YGrid','off')
xlim(ax2,[0 7])
hold off
    
    %% Cost breakdown

figure()
t = tiledlayout(1,1);
x = linspace(0.25,T-0.25,T*P);

bar(x,[output.results.distance_cost/1000 output.results.total_op_cost/1000], 'stacked')
legend('Fuel costs','Operation costs')
xlim([0 7])
ylim([0 200])
ax = gca;
ax.LineWidth = 0.5;
set(gca, 'YGrid', 'on', 'XGrid', 'on')
ax.XAxis.Visible = 'off';
ylabel('Costs [x10^3 MWK]')


ax2 = axes(t);
ax2.XAxisLocation = 'bottom';
ax2.YAxis.Visible  = 'off';
xticks(ax2, 0.5+[0 1 2 3 4 5 6 7])
xticklabels(ax2, {'Mon','Tue','Wed','Thu','Fri','Sat','Sun'});
ax2.Box = 'off';
set(ax2, 'XMinorGrid', 'off')

xlim(ax2,[0 7])


ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = 1;
grid minor
