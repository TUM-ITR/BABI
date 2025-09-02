function [MMH_samples, all_samples, acceptance_ratio] = staged_MMH(y, K, K_b, k_d, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, theta_init, log_pdf_x_0, x_0_init, proposal_cov_init, T_chunk, K_stage, alpha, regularizer)
%STAGED_MMH Run marginal Metropolis-Hastings (MMH) dsmplrt with incremental 
% data and adaptive proposal to obtain samples (theta, x_0:T) from the 
% joint parameter and state posterior distribution 
% p(theta, x_0:T | D=(y_training)).
% This method performs a staged Metropolis-Hastings (MH) procedure in which
% the number of data points used in the likelihood computation is gradually
% increased by a fixed chunk size. At each stage, the MMH sampler is run on 
% the current data subset, and the proposal distribution is adapted based 
% on the empirical covariance of the collected samples.
%
%   Inputs:
%       y: training output trajectory
%       K: number of models/scenarios to be sampled in final stage
%       K_b: length of the burn-in period applied at every stage
%       k_d: number of models/scenarios to be skipped to decrease correlation (thinning) in final stage
%       f_theta: state transition function parametrized by theta; has inputs (theta, x)
%       g_theta: measurement function parametrized by theta; has inputs (theta, x)
%       log_pdf_w_theta: function that returns the logarithm of the probability density function of the measurement noise parametrized by theta; has inputs (theta, w)
%       log_pdf_theta: function that returns the logarithm of the probability density function of theta (prior); has input (theta)
%       theta_init: initial theta
%       log_pdf_x_0: function that returns the logarithm of probability density function of x_0 (prior); has input (x_0)
%       x_0_init: initial x_0
%       proposal_cov_init: initial covariance of multivariate normal proposal distribution (e.g., covariance of the prior) 
%       T_chunk: number of datapoints added at each stage
%       K_stage: number of samples per stage
%       alpha: proposal scaling factor (scalar or vector containing it for each stage)
%       regularizer: small constant added to diagonal of proposal covariance
%
%   Outputs:
%       MMH_samples: final samples from full-data posterior
%       all_samples: cell array of samples per stage
%       acceptance_ratio: acceptance ratio of all stages
    
    % Get number of parameters, etc.
    n_theta = length(theta_init);
    n_x = length(x_0_init);
    n_z = n_theta + n_x;
    T = size(y, 2);

    N_stages = ceil(T / T_chunk);

    z_curr = [theta_init; x_0_init];
    proposal_cov = proposal_cov_init;
    all_samples = {};
    acceptance_ratio = zeros(1,N_stages);

    timer = tic;
    
    fprintf('### Started staged MMH sampling\n')
    for i = 1:N_stages
        % Select data chunk
        T_i = min(i * T_chunk, size(y, 2));
        y_i = y(:, 1:T_i);

        % Update proposal functions
        propose_z = @(z) mvnrnd(z', proposal_cov)';
        log_proposal_ratio_z = @(z, z_prop) 0;

        if i < N_stages
            % Run base MHH to generate K_stage samples without thinning
            [MMH_samples_stage, acceptance_ratio_stage] = MMH(y_i, K_stage, K_b, 0, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, z_curr(1:n_theta), log_pdf_x_0, z_curr(n_theta+1:end), propose_z, log_proposal_ratio_z, 'print', false);
        else
            % Run base MHH to generate K samples with thinning
            [MMH_samples_stage, acceptance_ratio_stage] = MMH(y_i, K, K_b, k_d, f_theta, g_theta, log_pdf_w_theta, log_pdf_theta, z_curr(1:n_theta), log_pdf_x_0, z_curr(n_theta+1:end), propose_z, log_proposal_ratio_z, 'print', false);
        end

        % Store samples
        all_samples{i} = MMH_samples_stage;
        acceptance_ratio(i) = acceptance_ratio_stage;
        
        if i < N_stages
            % Update proposal covariance
            Z = zeros(n_theta + n_x, K_stage);
            for j = 1:K_stage
                Z(:, j) = [MMH_samples_stage{j}.theta; MMH_samples_stage{j}.x_0];
            end
        
            cov_z = cov(Z');
            if isscalar(alpha)
                proposal_cov = alpha * cov_z + regularizer * eye(n_z);
            else
                proposal_cov = alpha(i) * cov_z + regularizer * eye(n_z);
            end

            z_curr = Z(:,end);
        else
            MMH_samples = MMH_samples_stage;
        end

        fprintf('Stage %i/%i complete\nAcceptance ratio: %.2f %%\n', i, N_stages, acceptance_ratio_stage)
    end
    average_acceptance_ratio = mean(acceptance_ratio);
    time_sampling = toc(timer);
    fprintf('### Staged MMH sampling complete\nRuntime: %.2f s\nAverage acceptance ratio: %.2f %%\n', time_sampling, average_acceptance_ratio);
end