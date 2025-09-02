function test_prediction(MMH_samples, n_x, f_theta, g_theta, sample_w_theta, k_n, y_test)
%TEST_PREDICTION Simulate the MMH samples forward in time and compare the predictions to the test data.
%
%   Inputs:
%       MMH_samples: MMH samples
%       n_x: number of states
%       f_theta: state transition function parametrized by theta; has inputs (theta, x, u)
%       g_theta: measurement function parametrized by theta; has inputs (theta, x, u)
%       sample_w_theta: function that returns a sample from the measurement noise distribution parametrized by theta; has input (theta)
%       k_n: each model is simulated k_n times
%       y_test: test output

fprintf('### Testing model\n')

% Get number of models, etc.
K = size(MMH_samples, 1);
n_y = size(y_test, 1);
T_test = size(y_test, 2);

% Pre-allocate.
x_test_sim = zeros(n_x, T_test+1, K, k_n);
y_test_sim = zeros(n_y, T_test, K, k_n);

%% Simulate models forward.
% If the Parallel Computing Toolbox is available, the for loop can be replaced with a parfor loop to decrease the runtime.
for k = 1:K
    % Get current model.
    f = @(x) f_theta(MMH_samples{k}.theta, x);
    g = @(x) g_theta(MMH_samples{k}.theta, x);
    sample_w = @() sample_w_theta(MMH_samples{k}.theta);

    % Simulate each model k_n times.
    for kn = 1:k_n
        % Pre-allocate.
        x_loop = zeros(n_x, T_test+1);
        y_loop = zeros(n_y, T_test);

        % Get initial state
        x_loop(:, 1) = MMH_samples{k}.x_t;

        % Simulate model forward.
        for t = 1:T_test
            x_loop(:, t+1) = f(x_loop(:, t));
            y_loop(:, t) = g(x_loop(:, t)) + sample_w();
        end

        % Store trajectory.
        x_test_sim(:, :, k, kn) = x_loop;
        y_test_sim(:, :, k, kn) = y_loop;
    end
end

% Reshape.
x_test_sim = reshape(x_test_sim, [n_x, T_test + 1, K * k_n]);
y_test_sim = reshape(y_test_sim, [n_y, T_test, K * k_n]);

fprintf('### Testing complete\n')

% Plot results.
plot_predictions(y_test_sim, y_test, 'plot_percentiles', true)

% Compute and print RMSE.
mean_rmse = sqrt(mean((squeeze(y_test_sim) - y_test').^2, 'all'));
fprintf('Mean rmse: %.2f\n', mean_rmse);
end