function [id, C, sumd] = kmeans_bisection(X, K, metric)
% Inputs: N x L data matrix X, where L is the number of data points and N
% is the number of variables in the system, scalar K representing the
% maximum number of clusters, and metric used to define distances (examples
% include 'sqeuclidean', 'cosine', and 'cityblock').
%
% Outputs: (K-1) x L matrix id, where the id(i,:) lists cluster labels of
% each observation with i+1 clusters chosen using the "bisecting" k-means
% algorithm; N x (2(K-1)) matrix C, where C(:,i) is the centroid position
% of the i+1th cluster; and 1 x (2(K-1)) array sumd, where sumd(i) is the
% combined point-to-centroid distance for cluster i.
%
% NOTE: For analysis in "Broken detailed balance and entropy production in
% the human brain" we use metric = 'cosine'

% Parameters:
maxIterations = 10000; % Maximum number of iterations for k-means algorithm
numRepeats = 1; % Number of repeats of k-means to perform at each level
split_criterion = 'spread'; % Pick next cluster to split based on spread

% Initialize variables:
N = size(X,1);
L = size(X,2);

id = zeros(K-1,L);
C = zeros(N, 2*(K-1));
sumd = zeros(1, 2*(K-1));

% Return if K == 1:
if K == 1
    return;
end

% Initialize temporary variables:
states_current = [0];
id_current = zeros(1,L);

% Loop over bisections of clusters:
for i = 2:K
    
    % Pick the next cluster to split:
    if length(states_current) == 1
        cluster = states_current;
    elseif strcmp(split_criterion, 'spread')
        [~, cInd] = max(sumd(states_current));
        cluster = states_current(cInd);
    elseif strcmp(split_criterion, 'size')
        cluster = mode(id_current);
    else
        return;
    end

    inds = find(id_current == cluster);
    X_current = X(:,inds);
    
    % Split chosen cluster using k-means:
    [id_temp, ~, sumd_temp] = kmeans(X_current', 2, 'Distance', metric,...
        'MaxIter', maxIterations, 'Replicates', numRepeats);
    id_temp = id_temp';
    sumd_temp = sumd_temp';
    
    % Define new cluster centroids:
    C_temp = zeros(N, 2);
    C_temp(:,1) = mean(X_current(:, id_temp == 1), 2);
    C_temp(:,2) = mean(X_current(:, id_temp == 2), 2);
    
    % Combine two new clusters with the old clusters:
    max_state = max(states_current);
    states_current = [states_current(states_current ~= cluster), (max_state + [1 2])];
    id_current(inds) = id_temp + max_state;
    
    % Save current variables:
    id(i-1,:) = id_current;
    sumd(max_state + [1 2]) = sumd_temp;
    C(:, max_state + [1 2]) = C_temp;
    
end


