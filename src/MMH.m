function [MMH_samples, acceptance_ratio] = MMH(y, K, K_b, k_d, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, theta_init, log_pdf_x_0, x_0_init, propose_z, log_proposal_ratio_z, varargin)
%MMH Run marginal Metropolis-Hastings (MMH) sampler to obtain samples (theta, x_0:T) from the joint parameter and state posterior distribution p(theta, x_0:T | D=(y_training)).
%
%   Inputs:
%       y: training output trajectory
%       K: number of models/scenarios to be sampled
%       K_b: length of the burn-in period (number of MH iterations)
%       k_d: number of MH iterations to be skipped between two stored samples to decrease correlation (thinning)
%       f_theta: state transition function parametrized by theta; has inputs (theta, x)
%       g_theta: measurement function parametrized by theta; has inputs (theta, x)
%       log_pdf_w_theta: function that returns the logarithm of the probability density function of the measurement noise parametrized by theta; has inputs (theta, w)
%       log_pdf_theta: function that returns the logarithm of the probability density function of theta (prior); has input (theta)
%       theta_init: initial theta
%       log_pdf_x_0: function that returns the logarithm of probability density function of x_0 (prior); has input (x_0)
%       x_0_init: initial x_0
%       propose_z: function that proposes new parameters z_prop = [theta_prop; x_0_prop] (proposal distribution); has input (z = [theta; x_0])
%       log_proposal_ratio_z: function that computes the logarithm of the ratio of proposal densities for z; has input arguments (z, z_prop)
%
%   Variable-length input argument list:
%       print: set to false to disable printout
%
%   Outputs:
%       MMH_samples: cell array containing the different MMH samples
%       acceptance_ratio: acceptance ratio

% Default values
print = true;

% Read variable-length input argument list.
for i = 1:2:length(varargin)
    if strcmp('print', varargin{i})
        print = varargin{i+1};
    end
end

%% Initialization
% Total number of MH iterations (length of the chain including the initial state)
N = K_b + 1 + (K - 1) * (k_d + 1);

% Get number of parameters, etc.
n_theta = length(theta_init);

% Time MMH sampling.
if print
    timer = tic;
end

% Initialize and pre-allocate
MMH_samples = cell(K, 1); % MMH samples
current_sample = 1;
n_accepted = 0;

% Initial state of the chain and corresponding log-likelihood and log-priors
theta = theta_init;
x_0 = x_0_init;
z = [theta; x_0];
log_p_theta = log_pdf_theta(theta);
log_p_x_0 = log_pdf_x_0(x_0);
log_likelihood = compute_log_likelihood(y, theta, x_0, f_theta, g_theta, log_pdf_w_theta);
log_target = log_likelihood + log_p_theta + log_p_x_0;

% Store initial state if no burn-in and thinning is used.
[MMH_samples, current_sample] = store_sample(MMH_samples, current_sample, 1, K_b, k_d, theta, x_0);

%% Run marginal Metropolis-Hastings sampler
for k = 1:N-1

    % Propose new parameters and a new initial state.
    z_prop = propose_z(z);
    theta_prop = z_prop(1:n_theta);
    x_0_prop = z_prop(n_theta+1:end);
    log_p_theta_prop = log_pdf_theta(theta_prop);
    log_p_x_0_prop = log_pdf_x_0(x_0_prop);

    % Proposals with zero prior density are rejected (the current state is repeated below).
    if isfinite(log_p_theta_prop) && isfinite(log_p_x_0_prop)
        % Forward simulation to determine likelihood of proposal.
        log_likelihood_prop = compute_log_likelihood(y, theta_prop, x_0_prop, f_theta, g_theta, log_pdf_w_theta);
        log_target_prop = log_likelihood_prop + log_p_theta_prop + log_p_x_0_prop;

        % Compute acceptance probability (always accept if the current state has zero or undefined density).
        if isfinite(log_target)
            log_acceptance_ratio = log_target_prop - log_target + log_proposal_ratio_z(z, z_prop);
        else
            log_acceptance_ratio = Inf;
        end

        % Accept or reject the proposal.
        if isfinite(log_target_prop) && rand(1) < exp(log_acceptance_ratio)
            n_accepted = n_accepted + 1;

            theta = theta_prop;
            x_0 = x_0_prop;
            z = z_prop;
            log_target = log_target_prop;
        end
    end

    % Sample k+1 of the chain is the current state (either the accepted proposal or the repeated previous state).
    % Use sample if the burn-in period is reached and the sample is not removed by thinning.
    [MMH_samples, current_sample] = store_sample(MMH_samples, current_sample, k+1, K_b, k_d, theta, x_0);

    % Print progress.
    if print && mod(k, max(1, floor((N - 1)/100))) == 0
        fprintf('%i/%i MH iterations completed\n', k, N-1)
    end
end

% Print runtime.
acceptance_ratio = (n_accepted / max(1, N - 1)) * 100;

if print
    time_sampling = toc(timer);
    fprintf('### MMH sampling complete\nRuntime: %.2f s\nAcceptance ratio: %.2f %%\n', time_sampling, acceptance_ratio)
end

end

function log_likelihood = compute_log_likelihood(y, theta, x_0, f_theta, g_theta, log_pdf_w_theta)
%COMPUTE_LOG_LIKELIHOOD Perform forward simulation starting at x_0 and compute the log-likelihood of the outputs y given theta and x_0.

T = size(y, 2);
x = x_0;
log_likelihood = 0;
for t = 1:T
    % Propagate state.
    if t >= 2
        x = f_theta(theta, x);
    end

    % Likelihood update based on measurement model (logarithms are used for numerical reasons).
    w = y(:, t) - g_theta(theta, x);
    log_likelihood = log_likelihood + log_pdf_w_theta(theta, w);
end

end

function [MMH_samples, current_sample] = store_sample(MMH_samples, current_sample, k, K_b, k_d, theta, x_0)
%STORE_SAMPLE Store sample k of the chain if it is past the burn-in period and not removed by thinning.

if (k > K_b) && (mod(k-(K_b + 1), k_d+1) == 0)
    MMH_samples{current_sample}.theta = theta;
    MMH_samples{current_sample}.x_0 = x_0;
    current_sample = current_sample + 1;
end

end
