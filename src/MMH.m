function [MMH_samples, acceptance_ratio] = MMH(y, K, K_b, k_d, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, theta_init, log_pdf_x_0, x_0_init, propose_z, log_proposal_ratio_z, varargin)
%MMH Run marginal Metropolis-Hastings (MMH) sampler to obtain samples (theta, x_0:T) from the joint parameter and state posterior distribution p(theta, x_0:T | D=(y_training)).
%
%   Inputs:
%       y: training output trajectory
%       K: number of models/scenarios to be sampled
%       K_b: length of the burn-in period
%       k_d: number of models/scenarios to be skipped to decrease correlation (thinning)
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
% Total number of samples to be generated
K_total = K_b + 1 + (K - 1) * (k_d + 1);

% Get number of parameters, etc.
n_theta = length(theta_init);
n_x = length(x_0_init);
T = size(y, 2);

% Time MMH sampling.
if print
    timer = tic;
end

% Initialize and pre-allocate
MMH_samples = cell(K, 1); % MMH samples
accepted_samples = 1;
current_sample = 1;
n_proposals = 0;
z = [theta_init; x_0_init];
log_likelihood = -Inf;
log_p_theta = -Inf;
log_p_x_0 = -Inf;

%% Run marginal Metropolis-Hastings sampler
while accepted_samples <= K_total

    % Propose new parameters and a new initial state.
    n_proposals = n_proposals + 1;
    z_prop = propose_z(z);
    theta_prop = z_prop(1:n_theta);
    x_0_prop = z_prop(n_theta+1:end);
    log_p_theta_prop = log_pdf_theta(theta_prop);
    log_p_x_0_prop = log_pdf_x_0(x_0_prop);
    if ~isfinite(log_p_theta_prop) || ~isfinite(log_p_x_0_prop)
        continue;
    end

    % Update model, i.e., update state transition and observation function and noise distribution.
    f = @(x) f_theta(theta_prop, x);
    g = @(x) g_theta(theta_prop, x);
    log_pdf_w = @(w) log_pdf_w_theta(theta_prop, w);

    x_prop = zeros(n_x, T);
    x_prop(:, 1) = x_0_prop;
    log_likelihood_prop = 0;

    % Forward simulation to determine likelihood of proposal.
    for t = 1:T
        % Propagate state.
        if t >= 2
            x_prop(:, t) = f(x_prop(:, t-1));
        end

        % Likelihood update based on measurement model (logarithms are used for numerical reasons).
        w = y(:, t) - g(x_prop(:, t));
        log_likelihood_prop = log_likelihood_prop + log_pdf_w(w);
    end

    % Compute acceptance probability.
    log_acceptance_ratio = log_likelihood_prop - log_likelihood ...
        +log_p_theta_prop - log_p_theta ...
        +log_p_x_0_prop - log_p_x_0 ...
        +log_proposal_ratio_z(z, z_prop);

    % Accept or reject the proposal.
    if rand(1) < exp(log_acceptance_ratio)
        accepted_samples = accepted_samples + 1;

        theta = theta_prop;
        x_0 = x_0_prop;
        z = [theta; x_0];
        log_likelihood = log_likelihood_prop;
        log_p_theta = log_p_theta_prop;
        log_p_x_0 = log_p_x_0_prop;

        % Use sample if the burn-in period is reached and the sample is not removed by thinning.
        if (K_b + 1 < accepted_samples) && (mod(accepted_samples-(K_b + 2), k_d+1) == 0)
            MMH_samples{current_sample}.theta = theta_prop;
            MMH_samples{current_sample}.x_0 = x_prop(:, 1);
            % MMH_samples{current_sample}.x_t = f(x_prop(:, end));
            current_sample = current_sample + 1;
        end

        % Print progress.
        if print
            fprintf('%i/%i samples accepted\n', accepted_samples-1, K_total)
        end
    end
end

% Print runtime.
acceptance_ratio = ((K_total - 1) / n_proposals) * 100;

if print
    time_sampling = toc(timer);
    fprintf('### MMH sampling complete\nRuntime: %.2f s\nAcceptance ratio: %.2f %%\n', time_sampling, acceptance_ratio)
end

end