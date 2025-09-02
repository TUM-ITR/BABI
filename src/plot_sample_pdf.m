function plot_sample_pdf(MMH_samples, varargin)
%PLOT_SAMPLE_PDF Plots the sample approximation of the pdf of the posterior via binning.
%
%   Inputs:
%       MMH_samples: MMH samples
%
%   Variable-length input argument list:
%       num_bins: number of bins

% Default values
num_bins = 100;

% Read variable-length input argument list.
for i = 1:2:length(varargin)
    if strcmp('num_bins', varargin{i})
        num_bins = varargin{i+1};
    end
end

% Get number of models.
K = size(MMH_samples, 1);

% Get parameters
n_theta = numel(MMH_samples{1}.theta);
theta_samples = zeros(n_theta, K);
for k = 1:K
    theta_samples(:, k) = MMH_samples{k}.theta;
end

% Get initial states
n_x = numel(MMH_samples{k}.x_0);
x_0_samples = zeros(n_x, K);
for k = 1:K
    x_0_samples(:, k) = MMH_samples{k}.x_0;
end

% Plot the pdf for all parameters.
for i = 1:n_theta
    figure();
    
    % Extract samples of this dimension
    s = theta_samples(i, :);
    
    % Create histogram and normalize to get PDF
    [counts, edges] = histcounts(s, num_bins, 'Normalization', 'pdf');
    
    % Bin centers for plotting
    bin_centers = 0.5 * (edges(1:end-1) + edges(2:end));
    
    % Plot
    bar(bin_centers, counts, 'FaceAlpha', 0.7);
    xlabel(['\theta_', num2str(i)]);
    ylabel('pdf');
    title(['p(\theta_', num2str(i), ' | D)']);
end

% Plot the pdf for all states.
for i = 1:n_x
    figure();
    
    % Extract samples of this dimension
    s = x_0_samples(i, :);
    
    % Create histogram and normalize to get PDF
    [counts, edges] = histcounts(s, num_bins, 'Normalization', 'pdf');
    
    % Bin centers for plotting
    bin_centers = 0.5 * (edges(1:end-1) + edges(2:end));
    
    % Plot
    bar(bin_centers, counts, 'FaceAlpha', 0.7);
    xlabel(['x_', num2str(i)]);
    ylabel('pdf');
    title(['p(x_', num2str(i), ' | D)']);
end
end