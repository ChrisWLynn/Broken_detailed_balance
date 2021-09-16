% Plot flux map. Use 'fluxMap_bootstrap' first. This script is specifically
% written for the rest HCP data

% Colormap for flux map:
color = .8*[181 23 0]/255;
X = linspace(1, 0, 1000)';
cmap = X*[1 1 1] + (1-X)*color;

% Increase scale of current:
scale = 4000; % For rest
% scale = 10000; % For rest (shuffle)
current_temp = scale*current;

% Units to measure current (in inverse seconds):
units = 5*10^(-4); % For rest
% units = 10^(-4); % For rest (shuffle)

% Number of standard deviations and alpha for ellipses:
num_stds = 2;
std_alpha = .2;

% Other plot options:
line_width = 1.5;
font_size = 21;
cmap_font_size = 18;

% Arrow appearance:
arrow_color = .1*[1 1 1];
arrow_width = line_width;
head_size = .1;

% Averages and standard deviations:
probability_mean = mean(probability,3);
divergence_mean = mean(current_divergence,3);
divergence_std = std(current_divergence,0,3);
current_mean = mean(current_temp,3);

% Calculate principle components and square roots of small and large
% eigenvalues of the covariance matrix for each current flux:
num_bins = size(current_mean,2);
current_std_position = position + current_mean;
current_std_comps = zeros(size(current_mean));
current_std_eigs = zeros(size(current_mean));

for i = 1:num_bins
    
    C = cov(reshape(current_temp(1,i,:), 1, num_bins), reshape(current_temp(2,i,:), 1, num_bins));
    [V, D] = eigs(C);
    [C_eigs, inds] = sort(diag(D));
    
    current_std_comps(:,i) = V(inds(1));
    current_std_eigs(:,i) = sqrt(C_eigs);
    
end 

% Make plot:

figure;
hold on;

% Plot histogram or divergence of flux vectors:
imagesc([bin_pos(1,2), bin_pos(1,end-1)], [bin_pos(2,2), bin_pos(2,end-1)], flipud(rot90(probability_mean)));

colormap(cmap);

c = colorbar;
c.Ticks = 0:.01:.1; % For flux map
lim = caxis;
% caxis([0 max(lim)])
caxis([0 .045]) % For comp with gambling
% caxis([0 .03]) % For shuffle
c.LineWidth = line_width;
c.FontSize = cmap_font_size;

% Add grid lines:
for i = 3:(size(bin_edges,2)-2)
    
    line(bin_edges(1,i)*[1 1], [bin_edges(2,1) bin_edges(2,end)], 'Color','w','LineStyle','-', 'LineWidth', 1);
    line([bin_edges(1,1) bin_edges(1,end)], bin_edges(2,i)*[1 1], 'Color','w','LineStyle','-', 'LineWidth', 1);
end

% Plot error ellipses for flux vectors:
for i = 1:num_bins
    
    angle = atan(-current_std_comps(2,i)/current_std_comps(1,i));
    plot_ellipse(num_stds*current_std_eigs(1,i), num_stds*current_std_eigs(2,i),...
        current_std_position(1,i), current_std_position(2,i),...
        angle, [0 0 0], std_alpha);
end

% Plot flux vectors:
quiver(position(1,:), position(2,:), current_mean(1,:), current_mean(2,:),...
    'Color', arrow_color, 'LineWidth', arrow_width, 'AutoScale', 'off', 'MaxHeadSize', head_size);

% Add key showing size of typical flux arrow:
fill([bin_pos(1, end-2) - .5*bin_size(1), bin_pos(1, end-1) + .5*bin_size(1), bin_pos(1, end-1) + .5*bin_size(1), bin_pos(1, end-2) - .5*bin_size(1)],...
    [bin_pos(2, end-2) - .5*bin_size(2), bin_pos(2, end-2) - .5*bin_size(2), bin_pos(2, end-1) + .5*bin_size(2), bin_pos(2, end-1) + .5*bin_size(2)],...
    'w', 'EdgeAlpha', 0);
quiver((bin_pos(1, end-2) - .3*bin_size(1))*[1 1], (bin_pos(2,end-2) - .3*bin_size(2))*[1 1],...
    [scale*units, 0], [0, scale*units], 'Color', arrow_color, 'LineWidth', arrow_width,...
    'AutoScale', 'off', 'MaxHeadSize', .5);
fill([bin_pos(1, end-2) - .5*bin_size(1), bin_pos(1, end-1) + .5*bin_size(1), bin_pos(1, end-1) + .5*bin_size(1), bin_pos(1, end-2) - .5*bin_size(1)],...
    [bin_pos(2, end-2) - .5*bin_size(2), bin_pos(2, end-2) - .5*bin_size(2), bin_pos(2, end-1) + .5*bin_size(2), bin_pos(2, end-1) + .5*bin_size(2)],...
    'w', 'EdgeAlpha', 0);
quiver((bin_pos(1, end-2) + .25*bin_size(1))*[1 1], (bin_pos(2,end-2) +.5*bin_size(2))*[1 1],...
    [scale*units, 0], [0, scale*units], 'Color', arrow_color, 'LineWidth', arrow_width,...
    'AutoScale', 'off', 'MaxHeadSize', .5);

xlabel('PC_1')
ylabel('PC_2')
ax = gca;
ax.XLim = [bin_edges(1,2), bin_edges(1,end-1)];
ax.YLim = [bin_edges(2,2), bin_edges(2,end-1)];
ax.XTick = -8:2:8;
ax.YTick = -8:2:8;
% ax.YTick = -500:250:500;
ax.TickLength = [0 0];
ax.LineWidth = line_width;
ax.FontSize = font_size;

set(gca, 'Layer', 'Top');
pbaspect([1 1 1])
box on;
hold off;