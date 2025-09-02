function plot_autocorrelation(MMH_samples, varargin)
%PLOT_AUTOCORRELATION Plot the autocorrelation function (ACF) of the MMH samples.
% This might be helpful when adjusting the thinning parameter k_d.
%
%   Inputs:
%       MMH_samples: MMH samples
%
%   Variable-length input argument list:
%       max_lag: maximum lag at which to calculate the ACF

% Get number of models.
K = size(MMH_samples, 1);

% Default values
max_lag = K - 1;

% Read variable-length input argument list.
for i = 1:2:length(varargin)
    if strcmp('max_lag', varargin{i})
        max_lag = varargin{i+1};
    end
end

% Get number of parameters of the MMH samples.
n_theta = numel(MMH_samples{1}.theta);
n_x = numel(MMH_samples{1}.x_0);

% Fill matrix with the series of the parameters of the MMH samples.
z_samples = zeros(K, n_theta + n_x);
for i = 1:K
    z_samples(i, :) = [MMH_samples{i}.theta; MMH_samples{i}.x_0];
end

% Calculate the autocorrelation.
autocorrelation = zeros(max_lag+1, n_theta + n_x);
for i = 1:(n_theta + n_x)
    [autocorrelation(:, i), ~] = autocorr(z_samples(:, i), 'NumLags', max_lag);
end

% Plot the ACF of the parameter samples.
figure();
hold on;
for i = 1:n_theta
        plot(0:max_lag, autocorrelation(:, i), 'LineWidth', 2, 'DisplayName', ['\theta_', num2str(i)]);
end

% Plot the ACF of the initial state samples.
% for i = 1:n_x        
%     plot(0:max_lag, autocorrelation(:, n_theta+i), 'LineWidth', 2, 'DisplayName', ['x_{0,', num2str(i),'}']);
% end

title('Autocorrelation Function (ACF)');
xlabel('Lag');
ylabel('ACF');
legend('show');
grid on;
ylim([-0.1, 1]);
yticks(-1:0.1:1);
xlim([0, 200]);
xticks(0:20:max_lag);
end