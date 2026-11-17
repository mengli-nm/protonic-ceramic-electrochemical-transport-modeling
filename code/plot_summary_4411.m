function plot_summary_4411()
% Plot the ensemble FE summaries for BZCYYb4411:
%   FE_All_4411.png : steam fractions at 50%, 70%, and 90%, four current densities.
%   RP_70_4411.png : 70% steam, three current densities, four Rp values.
% Input contains precomputed FE distributions for selected parameter ensembles.
% The input distributions are stored in figure_data/figure_summary_4411.mat.
base_dir = fileparts(mfilename('fullpath'));
input_mat = fullfile(base_dir, '..', 'figure_data', 'figure_summary_4411.mat');
output_dir = fullfile(base_dir, '..', 'generated');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
S = load(input_mat);

% FE vs current density at the three manuscript steam fractions.
fe_data = {S.FE_50, S.FE_70, S.FE_90};
fe_labels = {'50%', '70%', '90%'};
fe_colors = [20 32 238; 56 121 255; 133 203 245]/255;
current_fe = [-0.5, -1.0, -1.5, -2.0];
plot_grouped_boxes(fe_data, fe_labels, fe_colors, current_fe, [50 100], ...
    fullfile(output_dir, 'FE_All_4411.png'), 'Steam fraction');

% FE vs polarization resistance at 70% H2O only.
rp_data = {S.RP_70_0, S.RP_70_005, S.RP_70_01, S.RP_70_02};
rp_labels = {'R_p=0', 'R_p=0.05', 'R_p=0.1', 'R_p=0.2'};
rp_colors = [22 15 237; 55 87 240; 93 172 239; 155 251 250]/255;
current_rp = [-1.0, -1.5, -2.0];
plot_grouped_boxes(rp_data, rp_labels, rp_colors, current_rp, [40 100], ...
    fullfile(output_dir, 'RP_70_4411.png'), 'Polarization resistance');
end

function plot_grouped_boxes(group_data, group_names, group_colors, currents, ylim_values, output_file, legend_label)
% Draw boxes and whiskers directly so Statistics Toolbox is NOT needed.
nGroups = numel(group_data);
nCurr = numel(currents);
figure('Visible','off','Color','white'); hold on;
box_width = 0.16;
spacing = 0.21;
for k = 1:nGroups
    A = group_data{k}; % Each row is one current, entries are FE fractions.
    if size(A,1)~=nCurr, error('Incorrect FE summary dimensions.'); end
    for j = 1:nCurr
        xs = j + (k - (nGroups+1)/2)*spacing;
        draw_single_box(xs, 100*A(j,:), group_colors(k,:), box_width);
    end
end
set(gca,'XTick',1:nCurr,'XTickLabel',arrayfun(@(x)sprintf('%.1f',x),currents,'UniformOutput',false));
xlim([0.5, nCurr+0.5]); ylim(ylim_values);
xlabel('Current density (A cm^{-2})'); ylabel('Faradaic efficiency (FE, %)');
legend_handles = gobjects(nGroups,1);
for k=1:nGroups
    legend_handles(k)=patch(NaN,NaN,group_colors(k,:),'EdgeColor','black');
end
legend(legend_handles,group_names,'Location','southwest','Box','on');
set(gca,'FontSize',12,'Box','on');
saveas(gcf,output_file); close(gcf);
end

function draw_single_box(x, samples, color, box_width)
samples=sort(samples(isfinite(samples)));
if isempty(samples), return; end
q1=percentile_sorted(samples,0.25); med=percentile_sorted(samples,0.50); q3=percentile_sorted(samples,0.75);
iqr=q3-q1;
inliers=samples(samples >= q1-1.5*iqr & samples <= q3+1.5*iqr);
if isempty(inliers), inliers=samples; end
wlow=min(inliers); whigh=max(inliers);
patch([x-box_width/2 x+box_width/2 x+box_width/2 x-box_width/2], ...
      [q1 q1 q3 q3],color,'EdgeColor','black');
plot([x-box_width/2 x+box_width/2],[med med],'k-','LineWidth',1);
plot([x x],[wlow q1],'k-'); plot([x x],[q3 whigh],'k-');
plot([x-box_width/4 x+box_width/4],[wlow wlow],'k-');
plot([x-box_width/4 x+box_width/4],[whigh whigh],'k-');
outliers=samples(samples<wlow | samples>whigh);
if ~isempty(outliers), plot(x*ones(size(outliers)),outliers,'r+','MarkerSize',5); end
end

function q=percentile_sorted(a,p)
% Linear interpolation between ordered sample values.
n=numel(a); if n==1, q=a(1); return; end
idx=1+(n-1)*p; lo=floor(idx); hi=ceil(idx);
q=a(lo)+(a(hi)-a(lo))*(idx-lo);
end
