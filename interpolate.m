function X = interpolate(X0, X1, bin_edges)
% Inputs: Dx1 column vectors X0 and X1 defining beginning and ending points
% of the trajectory, where D is the dimensionality of the system, and
% dx(numBins + 1) matrix defining the edges of different bins discretizing
% the space.
%
% Output: DxL matrix X where the columns define unique positions in
% bin-space and X(:,i) and X(:,i+1) are directly adjacent.

% Dimensionality of system:
D = length(X0);

% Make linear interpolation between beginning and ending points:
steps = 0:.01:1;
X_space = X1*steps + X0*(1 - steps);

% Bin the interpolated trajectory:
X_bin = zeros(size(X_space));

for i = 1:D
    X_bin(i,:) = discretize(X_space(i,:), bin_edges(i,:));
end

% Loop through states to record transitions:
X = X_bin(:,1);

for i = 2:size(X_bin,2)
    
    X0_bin = X_bin(:,i-1);
    X1_bin = X_bin(:,i);
    
    if ~isequal(X1_bin, X0_bin)
        X = [X, X1_bin];
    end
end




