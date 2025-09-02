function delta = cp_upper_bound(V, K, beta)
%CP_UPPER_BOUND Compute the one-sided (1-beta) Clopper-Pearson upper confidence bound.
%
%   Inputs:
%     V: number of violations
%     K: number of samples in the test set
%     beta: confidence parameter
%
%   Output:
%     delta: upper bound in [0,1]

if V > K
    error('V must be <= K.');
end

if V == K
    delta = 1.0;
elseif V == 0
    % Zero violations
    delta = 1 - beta^(1 / K);
else
    % General case
    delta = betainv(1-beta, V+1, K-V);
end
end