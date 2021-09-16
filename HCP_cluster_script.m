% Script to perform hierarchical k-means clustering and estimate entropy
% productions for HCP data.

% List task scans:
tasks = {'EMOTION_LR'; 'EMOTION_RL'; 'GAMBLING_LR'; 'GAMBLING_RL'; 'LANGUAGE_LR';...
    'LANGUAGE_RL'; 'MOTOR_LR'; 'MOTOR_RL'; 'RELATIONAL_LR'; 'RELATIONAL_RL';...
    'REST1_LR'; 'REST1_RL'; 'SOCIAL_LR'; 'SOCIAL_RL'; 'WM_LR'; 'WM_RL'};

% Task groups (for which to estimate entropy production):
task_groups = {{'EMOTION_LR'; 'EMOTION_RL'}; {'GAMBLING_LR'; 'GAMBLING_RL'};...
    {'LANGUAGE_LR';'LANGUAGE_RL'}; {'MOTOR_LR'; 'MOTOR_RL'}; {'RELATIONAL_LR';...
    'RELATIONAL_RL'}; {'REST1_LR'; 'REST1_RL'}; {'SOCIAL_LR'; 'SOCIAL_RL'};...
    {'WM_LR'; 'WM_RL'}};

% Information to set:
inds_ID = 1:590; % Must be no greater than 590
num_clusters = 20; % Max number of clusters k
num_samples = 100; % Number of bootstrap samples

% Directory name:
dir_data = 'HCP_data_share';

% Information about tasks and subjects:
info = load([dir_data, '/information']);
num_regions = 100;

% Limit to specific set of tasks:
inds_task = find(ismember(info.task_info.task, tasks));
tasks = info.task_info.task(inds_task);
task_lengths = info.task_info.length(inds_task);

% Limit to specific set of subjects:
inds_complete = find(info.subject_info.complete);
IDs = info.subject_info.ID(inds_complete(inds_ID));

% Number of tasks, max number of clusters, and number of samples:
num_tasks = length(tasks);
num_groups = length(task_groups);
num_subjects = length(IDs);

% Create data structure:
for i = 1:num_tasks
    
    task = tasks{i};
    struct_temp = load([dir_data, '/HCP_data_task/', task]);
    
    % Add data and subject IDs:
    inds = ismember(struct_temp.task_IDs, IDs);
    data_struct.(task) = struct_temp.task_data(:,inds);
    ID_struct.(task) = struct_temp.task_IDs(inds);

end

% Keep the same number of data points for each task and each subject. Here
% we only consider the data at the beginning of tasks in order to treat
% each task equally:
task_length_min = min(task_lengths);
data_subsample = zeros(num_regions, num_tasks*task_length_min*num_subjects);
IDs_subsample = zeros(1, num_tasks*task_length_min*num_subjects);
tasks_subsample = cell(1, num_tasks*task_length_min*num_subjects);

for i = 1:num_tasks
    
    task = tasks{i};
    inds = ((i-1)*task_length_min*num_subjects + 1):(i*task_length_min*num_subjects);
    inds_beginning = subsample_data(ID_struct.(task), task_length_min);
    
    data_subsample(:, inds) = data_struct.(task)(:,inds_beginning);
    IDs_subsample(inds) = ID_struct.(task)(inds_beginning);
    tasks_subsample(inds) = repmat({task}, 1, length(inds));
    
end

% Data to save:
entProd = zeros(num_groups, num_clusters - 1, num_samples);
distance = zeros(num_clusters - 1, num_samples);
variance = zeros(num_clusters - 1, num_samples);
variance_tot = sum(var(data_subsample, [], 2));
centroids = cell(num_clusters - 1, num_samples);
missing_transitions = zeros(num_groups, num_clusters - 1, num_samples);
transition_mats = cell(num_groups, num_clusters - 1, num_samples);

% Loop over samples:
for i = 1:num_samples
    
    tic
    
    % Perform hierarchical clustering:
    [data_cluster, C, dist] = kmeans_bisection(data_subsample, num_clusters, 'cosine');
    
    % Loop over number of clusters:
    for j = 1:(num_clusters-1)
        
        data_cluster_temp = data_cluster(j,:);
        clusters = unique(data_cluster_temp);
        
        % Record cluster centroids:
        centroids{j,i} = C(:,clusters);
        
        % Calculate average distance from centroids:
        distance(j,i) = sum(dist(clusters))/length(data_cluster_temp);
        
        % Calculate variance:
        variance(j,i) = sum(var(C(:,data_cluster_temp), [], 2));
        
        % Bootstrap sample indices for transitions:
        samples = zeros(2, num_tasks*(task_length_min - 1)*num_subjects);
        tasks_sample = cell(1, num_tasks*(task_length_min - 1)*num_subjects);
        
        for k = 1:num_tasks
            
            inds_data = ((k-1)*task_length_min*num_subjects + 1):(k*task_length_min*num_subjects);
            inds_trans = ((k-1)*(task_length_min - 1)*num_subjects + 1):(k*(task_length_min - 1)*num_subjects);
            
            samples_task = bootstrap_transitions(inds_data, IDs_subsample(inds_data), 1, 1);

            samples(:, inds_trans) = samples_task;
            tasks_sample(inds_trans) = tasks_subsample(samples_task(1,:));
                
        end
        
        % Loop over task groups:
        for k = 1:num_groups
            
            % Get word samples for this task:
            inds_trans_task = find(ismember(tasks_sample, task_groups{k}));
            samples_task = samples(:, inds_trans_task);
            
            % Calculate entropy production and number of missing transitions
            [S, trans_missing, trans_mat] = entProd_transitions(data_cluster_temp, samples_task, 'naive');
            
            entProd(k,j,i) = S;
            missing_transitions(k,j,i) = trans_missing;
            transition_mats{k,j,i} = trans_mat;

        end
    end
    
    i
    toc
    
end

