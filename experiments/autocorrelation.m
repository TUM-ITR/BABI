% This script reproduces the results of the autocorrelation analysis
% reported in Section V of the paper (Figure 2).

% Clear.
clear;
clc;
close all;

% Import src.
addpath('..\src')

% Specify seed (for reproducible results).
rng(1);

%% Number of states, etc.
n_x = 2; % number of states
n_y = 1; % number of outputs

%% Learning parameters.
K_training = 1e5; % number of samples used to find barrier function
K_test = 0; % number of samples used to validate barrier function
K = K_training + K_test; % number of samples
k_d = 0; % number of samples to be skipped to decrease correlation (thinning)
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

%% Plot data.
% figure()
% hold on
% y_plot = plot(y, 'linewidth', 1);
% legend(y_plot, 'output')
% ylabel('y')
% xlabel('t')
%
% figure()
% hold on
% plot(x(1, :), 'linewidth', 1);
% plot(x(2, :), 'linewidth', 1);
% legend('x_1', 'x_2')
% ylabel('x')
% xlabel('t')

%% Learn models.
% Result: K models of the type
% x_t+1 = f_theta(MMH_samples{i}.theta, x_t, u_t)
% y_t = g_theta(MMH_samples{i}.theta, x_t, u_t) + w_t
% w_t ~ W_theta
[MMH_samples, ~, ~] = staged_MMH(y, K, K_b, k_d, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, theta_init, log_pdf_x_0, x_0_init, proposal_cov_init, T_chunk, K_stage, alpha, regularizer);

% Plot autocorrelation
plot_autocorrelation(MMH_samples, 'max_lag', 200);

% Plot sample pdf
% plot_sample_pdf(MMH_samples, 'num_bins', 500);
