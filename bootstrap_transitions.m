function samples = bootstrap_transitions(inds, IDs, order_sample, order_conserve)
% Inputs: Length L list of indices inds, ID assigned to each data point
% (i.e., each index in inds), order of words order_sample to sample from
% the data (e.g., order_sample = 1 returns simple transitions), and order
% of correlations order_conserve to conserve when sampling. For this
% function, we requre that order_conserve >= order_sample
%
% Output: Bootstrap sample words of length order_sample + 1 between the
% indices in inds. To preserve correlations of length order_conserve, we
% first sample words of length order_conserve + 1 and then pick out all
% words of length order_sample + 1. We return a list of words samples, were 
% samples(:,i) is the ith sample. The list is as long as number of words of
% length order_sample + 1 in the original list of indices inds. We only
% sample words corresponding to indices with the same ID.
%
% NOTE: For standard Markov analysis, we use order_sample = 1 and
% order_conserve = 1.

if order_sample > order_conserve
    error('Samples cannot be longer than correlations being conserved!');
end

L = length(inds);

% List of unique IDs:
IDs_unique = unique(IDs);
num_IDs = length(IDs_unique);

% Number of words:
num_words_sample = L - num_IDs*order_sample;
num_words_conserve = L - num_IDs*order_conserve;
samples_per_word = order_conserve - order_sample + 1;

% Make list of data points that are ok to sample:
inds_ok = zeros(1, num_words_conserve);
ind = 1;

for i = 1:num_IDs
    
    inds_ID = find(IDs == IDs_unique(i)); % For numerical IDs
%     inds_ID = find(strcmp(IDs, IDs_unique{i})); % For string IDs
    num_inds = length(inds_ID);
    
    inds_ok(ind:(ind + num_inds - order_conserve - 1)) = inds_ID(1:(end - order_conserve));
    
    ind = ind + num_inds - order_conserve;
    
end

% Sample beginnings of blocks of length order_conserve:
inds_start = randsample(inds_ok, ceil(num_words_sample/samples_per_word), true);

% Split each block into individual transitions:
samples = zeros(order_sample + 1, ceil(num_words_sample/samples_per_word)*samples_per_word);

count = 1;
for i = 1:length(inds_start)
    
    for j = 1:samples_per_word
        
        samples(:,count) = inds((inds_start(i) + j - 1):(inds_start(i) + j + order_sample - 1))';
        
        count = count + 1;
        
    end
end

samples = samples(:, 1:num_words_sample);

