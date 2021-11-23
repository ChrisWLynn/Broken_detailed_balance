function data_sample = bootstrap(data, IDs, order)
% Inputs: NxL data matrix data, ID assigned to each data point (i.e., each
% column of data), and order of samples to return order.
%
% Output: Bootstrap sampled data data_sample, where each sample has length
% order. Only sample consecutive data points corresponding to the same ID.

[N, L] = size(data);

% List of unique IDs:
IDs_unique = unique(IDs);
num_IDs = length(IDs_unique);

% Make list of data points that are ok to sample:
inds_ok = zeros(1, L - num_IDs*(order - 1));
ind = 1;

for i = 1:num_IDs
    
    inds_ID = find(IDs == IDs_unique(i));
    num_inds = length(inds_ID);
    
    inds_ok(ind:(ind + num_inds - order)) = inds_ID(1:(end - order + 1));
    
    ind = ind + num_inds - order + 1;
    
end

% Sample data points:
inds_sample = randsample(inds_ok, ceil(L/order), true);

data_sample = zeros(N, order*ceil(L/order));

for i = 1:length(inds_sample)
    
    data_sample(:, (order*(i-1) + 1):(order*i)) = data(:, inds_sample(i):(inds_sample(i) + order - 1));
    
end

data_sample = data_sample(:,1:L);