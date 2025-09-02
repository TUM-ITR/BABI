% This helper function runs a simulation of the baseline method (sampling
% from the prior instead of the posterior) for the passed seed.

function result = run_simulation_baseline(seed)
% Specify seed (for reproducible results).
rng(seed);

%% Number of states.
n_x = 2; % number of states

%% Safety certificate parameters.
beta = 0.001; % confidence parameter

gamma = .01;
r_min = 1;
r_max = 3;
order_h = [0, 2, 4, 6]; % polynomial order for barrier function
order_q123 = [0, 2, 4]; % polynomial order for functions q1,q2,q3
order_q45 = 2; % polynomial order for functions q4,5
orders = {order_h, order_q123, order_q45};

% Bounds on acceptable feasibility ratio. Runs outside this interval are counted as failed.
min_feasability_ratio = 0.9; % lower bound
max_feasability_ratio = 1.1; % upper bound

%% Learning parameters.
K_training = 750; % number of samples used to find barrier function
K_test = 750; % number of samples used to validate barrier function
K = K_training + K_test; % number of samples

%% State-space prior and proposal distributions.
% Define basis functions to simplify notation
phi = @(x) [x(1); x(2); x(1)^3; 1];

% State transition function
dt = .1;
f_theta = @(theta, x) x + dt * [1, -1, -1, 0; theta(1), theta(2), 0, theta(3)] * phi(x);

% Normally distributed prior for parameters
theta_mean = [0.0; 0.0; 0.0];
theta_var = [.1; .1; .1];
theta_cov = diag(theta_var);

% Normally distributed prior for x_0
x_0_mean = [0; 0];
x_0_var = [.1; .1];
x_0_cov = diag(x_0_var);

%% Sample parameters from the prior.
prior_samples = cell(K, 1);

for i = 1:K
    prior_samples{i}.theta = mvnrnd(theta_mean, theta_cov)';
    prior_samples{i}.x_0 = mvnrnd(x_0_mean, x_0_cov)';
end

% Plot autocorrelation
% plot_autocorrelation(prior_samples, 'max_lag', 200);

% Plot sample pdf
% plot_sample_pdf(prior_samples, 'num_bins', 100);

% Split into training and validation set
prior_samples_training = prior_samples(1:K_training);
prior_samples_test = prior_samples(K_training+1:end);

%% Find barrier function candidate.
% Solve problem with all samples.
fprintf("### Started search for barrier function\n");
[h, q1, q2, q3, q4, q5, info] = find_candidate(prior_samples_training, f_theta, n_x, r_min, r_max, gamma, orders);

%% Validate the barrier function candidate on a test set and check if it is valid for the true system.
if (info.numerr == 0) && (info.pinf == 0) && (info.dinf == 0) && (info.numerr == 0) && (info.feasratio >= min_feasability_ratio) && (info.feasratio <= max_feasability_ratio)
    fprintf("### Barrier function candidate found. Starting validation.\n");
    [delta, ~] = validate_candidate(beta, prior_samples_test, f_theta, n_x, r_min, r_max, gamma, h, q1, q2, q3, q4, q5, 'min_feasability_ratio', min_feasability_ratio, 'max_feasability_ratio', max_feasability_ratio);

    % Check if the candidate is valid for the true system.
    true_system = cell(1,1);
    true_system{1}.theta = theta_true;
    true_system{1}.x_0 = x(:, 1);

    [~, V] = validate_candidate(0, true_system, f_theta, n_x, r_min, r_max, gamma, h, q1, q2, q3, q4, q5, 'min_feasability_ratio', min_feasability_ratio, 'max_feasability_ratio', max_feasability_ratio);

    % Candidate is valid if there is no violation for the true system.
    if V == 0
        valid = 1;
    else
        valid = 0;
    end
else
    delta = NaN;
    valid = NaN;
end

%% Return the results.
result.seed = seed;
result.prior_samples = prior_samples;
result.h = h;
result.q1 = q1;
result.q2 = q2;
result.q3 = q3;
result.q4 = q4;
result.q5 = q5;
result.info = info;
result.delta = delta;
result.valid = valid;
end