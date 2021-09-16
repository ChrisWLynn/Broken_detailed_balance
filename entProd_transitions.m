function [S, missing_transitions, C] = entProd_transitions(X, inds, type)
% Input: NxL data array X where each column X(:,i) is a different
% observation and 2x(L-1) vector inds where ind(:,i) are the indices of
% the i^th transition to consider. We also take in a string type, which can
% either be 'naive', 'add1', or 'add2' corresponding to no correction,
% adding 1 to each entry of the backward counts, and adding one entry to
% both forward and backward counts.
%
% Output: Estimate of entropy production S for Markov processes (that is,
% only considering words of length 2). Also return the number of missing
% transitions 'missing_transitions' and the transition counts matrix C.
%
% NOTE: For the analysis in "Broken detailed balance and entropy production
% in the human brain" we use type = 'naive'

if ~ismember(type, {'naive', 'add1', 'add2'})
    error('type must take value naive, add1, or add2!');
end

% Turn state observations into numbers:
[~,~,X_new] = unique(X', 'rows');
N = max(X_new);

% Calculate entropy production:
 
C = zeros(N);

for i = 1:size(inds,2)
    
    C(X_new(inds(1,i)), X_new(inds(2,i))) = C(X_new(inds(1,i)), X_new(inds(2,i))) + 1;
    
end

missing_transitions = sum(C(:) == 0);

% Calculate forward and backward probabilities:
if strcmp(type, 'naive')
    
    P_f = C/sum(C(:));
    P_b = P_f';
    
elseif strcmp(type, 'add1')
    
    P_f = C/sum(C(:));
    P_b = (C + 1)'/sum(sum(C + 1));
    
elseif strcmp(type, 'add2')
    
    P_f = (C + 1)/sum(sum(C + 1));
    P_b = P_f';
    
end

% Calculate entropy production:
logP = log2(P_f./P_b);
logP(P_f == 0) = 0;
logP(P_b == 0) = 0;

S = sum(sum(P_f.*logP));
