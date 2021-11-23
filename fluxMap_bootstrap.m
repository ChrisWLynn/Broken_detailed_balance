% Script to calculate flux map of brain activity from HCP data.

% Parameters to set:
num_samples = 100; % Number of bootstrap samples
order = 1; % Order of correlations preserved in bootstrap samples (keep 1 for Markov transitions)
num_comps = 2; % Number of principle components
num_bins = 10; % Number of discrete bins to consider in each dimension
num_STD = 2; % Number of standard deviations from component center to cut off flux analysis

% Length of time step in seconds:
dt = .72;

% Directory:
directory = 'HCP_data_share/HCP_data_task/';

% HCP data:
struct1_LR = load([directory, 'REST1_LR']);
struct1_RL = load([directory, 'REST1_RL']);
struct2_LR = load([directory, 'GAMBLING_LR']);
struct2_RL = load([directory, 'GAMBLING_RL']);

X_data = [struct1_LR.task_data, struct1_RL.task_data, struct2_LR.task_data, struct2_RL.task_data];
ID_data = [struct1_LR.task_IDs, struct1_RL.task_IDs + max(struct1_LR.task_IDs),...
    struct2_LR.task_IDs + max(struct1_LR.task_IDs) + max(struct1_RL.task_IDs),...
    struct2_RL.task_IDs + max(struct1_LR.task_IDs) + max(struct1_RL.task_IDs) + max(struct2_LR.task_IDs)];
task_data = [ones(1, length(struct1_LR.task_IDs) + length(struct1_LR.task_IDs)),...
    2*ones(1, length(struct2_LR.task_IDs) + length(struct2_LR.task_IDs))];

% List of different IDs and tasks/settings:
ID_unique = unique(ID_data);
num_IDs = length(ID_unique);

% Restrict to a certain amount of data per subject:
L_ID = 176;

X_data_temp = zeros(size(X_data,1), L_ID*num_IDs);
ID_data_temp = zeros(1, L_ID*num_IDs);
task_data_temp = zeros(1, L_ID*num_IDs);

for i = 1:num_IDs
    
    inds = find(ID_data == ID_unique(i));
    X_data_temp(:, (L_ID*(i-1) + 1):(L_ID*i)) = X_data(:, inds(1:L_ID));
    ID_data_temp((L_ID*(i-1) + 1):(L_ID*i)) = ID_data(inds(1:L_ID));
    task_data_temp((L_ID*(i-1) + 1):(L_ID*i)) = task_data(inds(1:L_ID));
    
end

X_data = X_data_temp;
ID_data = ID_data_temp;
task_data = task_data_temp;

% Total length of data:
t_tot = dt*size(X_data,2);

% Z-score data:
X_data = (X_data - mean(X_data(:)))/std(X_data(:));

% Perform PCA and bin data:
[PCs, X_PCA, ~] = pca(X_data', 'NumComponents', num_comps);
X_PCA = X_PCA';

% Bin PCA data:
X_bin = zeros(size(X_PCA));
bin_pos = zeros(num_comps, num_bins + 2);
bin_size = zeros(num_comps,1);
bin_edges = zeros(num_comps, num_bins + 3);

for i = 1:num_comps
    
    mean_comp = mean(X_PCA(i,:));
    std_comp = std(X_PCA(i,:));
    
    bin_size(i) = 2*std_comp*num_STD/num_bins;
    
    bin_edges(i,:) = linspace(mean_comp - num_STD*std_comp - bin_size(i),...
        mean_comp + num_STD*std_comp + bin_size(i), num_bins + 3);
    
    bin_pos(i,:) = (bin_edges(i,1:(end-1)) + bin_edges(i,2:end))/2;
    
    X_bin(i,:) = discretize(X_PCA(i,:), bin_edges(i,:));
end

%% Restrict to given task/setting:

% Task/setting for which to calculate flux currents (here, 1 for rest or 2
% for gambling task):
task_choice = 1; 

% Only keep data for given task/setting:
inds_task = find(task_data == task_choice);
X_PCA_task = X_PCA(:,inds_task);
X_bin_task = X_bin(:,inds_task);
ID_data_task = ID_data(inds_task);
ID_unique = unique(ID_data_task);

% Create list of state transitions:
transitions = zeros(5, num_bins*length(inds_task));

count = 1;

for i = 1:length(ID_unique)
    
    inds_ID = find(ID_data_task == ID_unique(i));
    X_bin_temp = X_bin_task(:,inds_ID);
    X_PCA_temp = X_PCA_task(:,inds_ID);

    for j = 2:size(X_bin_temp,2)
        
        X0 = X_bin_temp(:,j-1);
        X1 = X_bin_temp(:,j);
        
        % Test if X1 falls within flux map:
        if ~isnan(sum(X1))
            
            % Test if X0 falls within flux map:
            if ~isnan(sum(X0))
                
                % Test if X1 is directly adjacent to X0:
                if (ismember(X1(1), X0(1) + [-1 0 1]) && (X1(2) == X0(2))) ||...
                        ((X1(1) == X0(1)) && ismember(X1(2), X0(2) + [-1 0 1]))
                    
                    transitions(:, count) = [X0; X1; ID_unique(i)];
                    count = count + 1;
                    
                else
                    
                    % Linear interpolation between non-adjacent states:
                    X_temp = interpolate(X_PCA_temp(:,j-1), X_PCA_temp(:,j), bin_edges);
                    
                    for k = 2:size(X_temp,2)
                        X0_temp = X_temp(:,k-1);
                        X1_temp = X_temp(:,k);
                        
                        transitions(:, count) = [X0; X1; ID_unique(i)];
                        count = count + 1;
                        
                    end
                    
                end             
            end
        end
        
    end
    
end

transitions = transitions(:, transitions(5,:) > 0);

% Variables to keep track of:
mapSize = num_bins*ones(1, num_comps);
probability = zeros([mapSize, num_samples]); % State probabilities in matrix format
current_divergence = zeros([mapSize, num_samples]); % Divergence of probability current in matrix format
position = zeros(num_comps, num_bins^num_comps); % Bin positions in quiver format
current = zeros(num_comps, num_bins^num_comps, num_samples); % Probability current in quiver format

% Record bin positions in quiver format:
ind = 0;
for i = 2:(num_bins+1)
    for j = 2:(num_bins+1)
        
        ind = ind + 1;
        position(:,ind) = [bin_pos(1,i); bin_pos(2,j)];
        
    end
end
        
% Loop over bootstrap samples:
for i = 1:num_samples
    
    transitions_sample = bootstrap(transitions(1:4,:), transitions(5,:), order);
    
    % Loop through data and count state transitions:
    mapSize_temp = (num_bins+2)*ones(1, num_comps);
    prob_sample = zeros(mapSize_temp);
    trans_rate_sample = zeros([mapSize_temp, mapSize_temp]);
    
    % Count first state in probability:
    prob_sample(transitions_sample(1:2, 1)) = 1;
    
    for j = 1:size(transitions_sample,2)
        
        % Count new state in probability:
        prob_sample(transitions_sample(3,j), transitions_sample(4,j)) =...
            prob_sample(transitions_sample(3,j), transitions_sample(4,j)) + 1;
        
        % Count transition:
        trans_rate_sample(transitions_sample(1,j), transitions_sample(2,j),...
            transitions_sample(3,j), transitions_sample(4,j)) =...
            trans_rate_sample(transitions_sample(1,j), transitions_sample(2,j),...
            transitions_sample(3,j), transitions_sample(4,j)) + 1;
    
    end
    
    % Normalize:
    probability(:,:,i) = prob_sample(2:(end-1), 2:(end-1))/sum(prob_sample(:));
    trans_rate_sample = trans_rate_sample/t_tot;
    
    % Calculate probability flux:
    trans_mat = reshape(trans_rate_sample, [(num_bins+2)^num_comps, (num_bins+2)^num_comps]);
    flux_rate = reshape(trans_mat - trans_mat', size(trans_rate_sample));
    
    % Convert probability flux into quiver format and calculate flux divergence:
    ind = 0;
    for j = 2:(num_bins+1)
        for k = 2:(num_bins+1)
            
            ind = ind + 1;
            
            % Calculate probability flux vector:
            current(:, ind, i) = 1/2*[flux_rate(j-1, k, j, k) + flux_rate(j, k, j+1, k);...
                flux_rate(j, k-1, j, k) + flux_rate(j, k, j, k+1)];
            
            % Calculate divergence of probability flux:
            current_divergence(j-1,k-1,i) = flux_rate(j-1, k, j, k) - flux_rate(j, k, j+1, k) +...
                flux_rate(j, k-1, j, k) - flux_rate(j, k, j, k+1);
            
        end
    end
    
end
