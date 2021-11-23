function inds = subsample_data(labels, N)
% Inputs: Array of labels and number of indices N to pick for each label
%
% Output: For each label, choose N consecutive indices. Return the total
% list of indices inds.

labels_unique = unique(labels);
num_labels = length(labels_unique);

inds = zeros(1, num_labels*N);

% Loop over different labels:
for i = 1:num_labels
    
    inds_label = find(labels == labels_unique(i));
    
    % If there are not enough indices, return nothing:
    if length(inds_label) < N
        inds = [];
        return;
    end
    
    % Otherwise, subsample the indices:
%     ind_temp = randi(length(inds_label) - N + 1); % Start at random point
    ind_temp = 1; % Start at beginning
%     ind_temp = length(inds_label) - N + 1; % Start at the end
    
    inds((N*(i-1) + 1):(N*i)) = inds_label(ind_temp:(ind_temp + N - 1));
    
end