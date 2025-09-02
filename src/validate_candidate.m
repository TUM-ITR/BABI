function [delta, V] = validate_candidate(beta, MMH_samples_test, f_theta, n_x, r_min, r_max, gamma, h, q1, q2, q3, q4, q5, varargin)
%VALIDATE Test barrier function candidate on an independent test set and
% compute bound on the probability that the candiate fails to be forward
% invariant for the true, unknown system.
%
%   Inputs:
%       beta: confidence parameter
%       MMH_samples_test: samples in the test set
%       f_theta: state transition function parametrized by theta; has inputs (theta, x)
%       n_x: number of states
%       r_min: radius of the circular set X_min centered at the origin
%       r_max: radius of the circular set X_max centered at the origin
%       gamma: parameter in the class K function alpha(x) = gamma * x
%       h: candidate barrier function
%       q1: function q_1
%       q2: function q_2
%       q3: function q_3
%       q4: function q_4
%       q5: function q_5
%
%   Variable-length input argument list:
%       min_feasability_ratio: lower bound on acceptable feasibility ratio. Runs below this threshold are counted as failed.
%       max_feasability_ratio: upper bound on acceptable feasibility ratio. Runs above this threshold are counted as failed.
%
%   Outputs:
%       delta: upper bound on the probability that h fails to be forward invariant
%	    V: number of violations

% Default values
min_feasability_ratio = 0.9;
max_feasability_ratio = 1.1;

% Read variable-length input argument list.
for i = 1:2:length(varargin)
    if strcmp('min_feasability_ratio', varargin{i})
        min_feasability_ratio = varargin{i+1};
    elseif strcmp('max_feasability_ratio', varargin{i})
        max_feasability_ratio = varargin{i+1}; 
    end
end

%% Initialize.
K_test = numel(MMH_samples_test);
V = 0;

for k = 1:K_test
    %% Set up SOS program
    x = mpvar('x', n_x, 1);
    prog = sosprogram(x);

    %% Add constraints
    relaxation1 = r_max^2 - (x' * x);
    f = f_theta(MMH_samples_test{k}.theta, x);
    expr = subs(h, x, f) - (1 - gamma) * h;
    prog = sosineq(prog, expr-q1*relaxation1);

    h_min = subs(q4, x(1), x'*x) - subs(q4, x(1), r_min^2);
    h_max = subs(q5, x(1), x'*x) - subs(q5, x(1), r_max^2);
    expr2 = h_max - h;
    relaxation2 = r_min^2 - (x' * x);
    expr3 = h - h_min - .001;

    prog = sosineq(prog, expr2-q2*relaxation1);
    prog = sosineq(prog, expr3-q3*relaxation2);

    %% Solve SOS program
    solver_opt.solver = 'mosek';
    solver_opt.simplify = 'on';

    % solver_opt.params.MSK_IPAR_NUM_THREADS = 1;
    % solver_opt.params.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
    % solver_opt.params.MSK_DPAR_INTPNT_TOL_MU_RED = 1e-10;

    prog = sossolve(prog, solver_opt);

    %% Check if SOS program is feasible
    if (prog.solinfo.info.numerr ~= 0) || (prog.solinfo.info.pinf ~= 0) || (prog.solinfo.info.dinf ~= 0) || (prog.solinfo.info.feasratio < min_feasability_ratio) || (prog.solinfo.info.feasratio > max_feasability_ratio)
        V = V + 1;
    end
end

%% Compute bound delta
delta = cp_upper_bound(V, K_test, beta);
end