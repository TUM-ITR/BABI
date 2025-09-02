% This script reproduces the results of the exemplary run reported in
% Section V of the paper (Figure 3).

% Clear.
clear;
clc;
close all;

% Import src.
addpath('..\src')

% Import SOSTOOLS and Mosek.
addpath(genpath('Path\to\SOSTOOLS'));
addpath(genpath('Path\to\MOSEK'));

% Specify seed (for reproducible results).
rng(1);

%% Number of states, etc.
n_x = 2; % number of states
n_y = 1; % number of outputs

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
k_d = 150; % number of samples to be skipped to decrease correlation (thinning)
K_b = 200; % length of burn-in period for each stage
T_chunk = 10; % number of datapoints added at each stage
K_stage = 500; % number of samples per stage
alpha = linspace(1, 0.25, 200);
regularizer = 0; % regularizer for proposal covariance

%% State-space prior and proposal distributions.
% Define basis functions to simplify notation
phi = @(x) [x(1); x(2); x(1)^3; 1];

% State transition function
n_theta = 3; % number of parameters
dt = .1;
f_theta = @(theta, x) x + dt * [1, -1, -1, 0; theta(1), theta(2), 0, theta(3)] * phi(x);

% Measurement model - assumed to be known (not dependent on theta)
g_theta = @(theta, x) [1, 0] * x;
% Zero-mean Gaussian measurement noise with variance R - normalizing factors are ommited as they cancel out in the acceptance ratio
R = 0.1;
log_pdf_w_theta = @(theta, w) -(w.^2) / (2 * R);
sample_w_theta = @(theta) normrnd(0, sqrt(R));

% Normally distributed prior for parameters
theta_mean = [0.0; 0.0; 0.0];
theta_var = [.1; .1; .1];
theta_cov = diag(theta_var);
% Log pdf - normalizing factors are ommited as they cancel out in the acceptance ratio
log_pdf_theta = @(theta) -0.5 * (theta - theta_mean)' * (theta_cov \ (theta - theta_mean));
theta_init = theta_mean;

% Normally distributed prior for x_0
x_0_mean = [0; 0];
x_0_var = [.1; .1];
x_0_cov = diag(x_0_var);
% Log pdf - normalizing factors are ommited as they cancel out in the acceptance ratio
log_pdf_x_0 = @(x_0) -0.5 * (x_0 - x_0_mean)' * (x_0_cov \ (x_0 - x_0_mean));
x_0_init = x_0_mean;

proposal_cov_init = diag([theta_var; x_0_var]);

%% Parameters for data generation.
T = 2000; % number of steps for training

%% Generate training data.
% Choose the actual system (to be learned) and generate input-output data f length T_all.
% The system is of the form
% x_t+1 = f_true(x_t)
% y_t = g_true(x_t) + N(0, R_true).

% True, unknown system.
theta_true = [0.01;-0.005;0];

f_true = @(x) f_theta(theta_true, x); % true state transition function
g_true = @(x) g_theta(theta_true, x); % true measurement function
R_true = R; % true variance of zero-mean Gaussian measurement noise

% Generate data by forward simulation.
x = zeros(n_x, T+1); % true latent state trajectory
x(:, 1) = mvnrnd(x_0_mean', diag(x_0_var))'; % random initial state
y = zeros(n_y, T); % output trajectory (measured)

for t = 1:T
    x(:, t+1) = f_true(x(:, t));
    y(:, t) = g_true(x(:, t)) + mvnrnd(zeros(n_y,1), R_true)';
end

%% Learn models.
% Result: K models of the type
% x_t+1 = f_theta(MMH_samples{i}.theta, x_t, u_t)
% y_t = g_theta(MMH_samples{i}.theta, x_t, u_t) + w_t
% w_t ~ W_theta
[MMH_samples, ~, ~] = staged_MMH(y, K, K_b, k_d, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, theta_init, log_pdf_x_0, x_0_init, proposal_cov_init, T_chunk, K_stage, alpha, regularizer);

% Plot autocorrelation
% plot_autocorrelation(MMH_samples, 'max_lag', 200);

% Plot sample pdf
% plot_sample_pdf(MMH_samples, 'num_bins', 100);

% Split into training and validation set
MMH_samples_training = MMH_samples(1:K_training);
MMH_samples_test = MMH_samples(K_training+1:end);

%% Find barrier function candidate.
% Solve problem with all samples.
fprintf("### Started search for barrier function\n");
[h, q1, q2, q3, q4, q5, info] = find_candidate(MMH_samples_training, f_theta, n_x, r_min, r_max, gamma, orders);

%% Validate the barrier function candidate on a test set and check if it is valid for the true system.
if (info.numerr == 0) && (info.pinf == 0) && (info.dinf == 0) && (info.numerr == 0) && (info.feasratio >= min_feasability_ratio) && (info.feasratio <= max_feasability_ratio)
    fprintf("### Barrier function candidate found. Starting validation.\n");
    [delta, ~] = validate_candidate(beta, MMH_samples_test, f_theta, n_x, r_min, r_max, gamma, h, q1, q2, q3, q4, q5, 'min_feasability_ratio', min_feasability_ratio, 'max_feasability_ratio', max_feasability_ratio);

    % Check if the candidate is valid for the true system.
    true_system = cell(1,1);
    true_system{1}.theta = theta_true;
    true_system{1}.x_0 = x(:, 1);

    [~, V] = validate_candidate(0, true_system, f_theta, n_x, r_min, r_max, gamma, h, q1, q2, q3, q4, q5, 'min_feasability_ratio', min_feasability_ratio, 'max_feasability_ratio', max_feasability_ratio);

    % Candidate is valid if there is no violation for the true system.
    if V == 0
        valid = "yes";
    else
        valid = "no";
    end
    feasible = "yes";

    % Plot the forward invariant set.
    % Simulate the true system for a long time such the trajectory
    % approaches the limit cycle.
    x_limit_cycle = zeros(n_x, 10*T+1);
    x_limit_cycle(:, 1) = mvnrnd(x_0_mean', diag(x_0_var))';
    for t = 1:10*T
        x_limit_cycle(:, t+1) = f_true(x_limit_cycle(:, t));
    end
    x_limit_cycle = x_limit_cycle(:,5*T:end);
    plot_invariant_set(h, f_true, r_min, r_max, 'x_trajectory', x_limit_cycle, 'x_trajectory_label', 'Limit cycle');
else
    delta = NaN;
    valid = "N/A";
    feasible = "no";
end

fprintf("### Simulation complete.\n")
fprintf("SOS program feasible: %s\n", feasible)
fprintf("1-delta: %.2f %%\n", 100*(1-delta))
fprintf("Valid for true system: %s\n", valid)
fprintf("Barrier function h(x):\n")
h