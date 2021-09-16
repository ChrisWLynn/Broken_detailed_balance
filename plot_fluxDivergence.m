% Plot divergence of probability current. Use 'fluxMap_bootstrap' first.

% Colormap:
X = linspace(1, 0, 1000)';
cmap = [(1-X)*[1 1 1] + X*[0 77 128]/255; X*[1 1 1] + (1-X)*.9*[150 23 0]/255];

% Averages and standard deviations:
divergence_mean = mean(current_divergence,3);
divergence_std = std(current_divergence,0,3);

% Make plot:

figure;
hold on;

% Plot divergence of flux vectors:
imagesc([bin_pos(1,2), bin_pos(1,end-1)], [bin_pos(2,2), bin_pos(2,end-1)], flipud(rot90(divergence_mean./divergence_std)));

colormap(cmap);

c = colorbar;
c.Ticks = -.8:.2:.8;
c.LineWidth = 1;
lim = caxis;
caxis([-.85 .85])
c.LineWidth = 1.5;

% Add grid lines:
for i = 3:(size(bin_edges,2)-2)
    
    line(bin_edges(1,i)*[1 1], [bin_edges(2,1) bin_edges(2,end)], 'Color','w','LineStyle','-', 'LineWidth', 1);
    line([bin_edges(1,1) bin_edges(1,end)], bin_edges(2,i)*[1 1], 'Color','w','LineStyle','-', 'LineWidth', 1);
end

xlabel('PC_1')
ylabel('PC_2')
ax = gca;
ax.XLim = [bin_edges(1,2), bin_edges(1,end-1)];
ax.YLim = [bin_edges(2,2), bin_edges(2,end-1)];
% ax.YTick = -500:250:500;
ax.TickLength = [0 0];
ax.LineWidth = 1.5;
ax.FontSize = 21;

set(gca, 'Layer', 'Top');
pbaspect([1 1 1])
box on;
hold off;