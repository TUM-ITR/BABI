% This script reproduces the results of the Monte Carlo simulation 
% of the propsed method reported in Section V of the paper 
% (Table 1; column "Posterior (proposed method)"). 
% It calls the function run_simulation with seeds 1:100 and stores the 
% results. The Parallel Computing Toolbox should be used to parallelize 
% the simulations.

% Clear.
clear;
clc;
close all;

% Import src.
addpath('..\src')

% Import SOSTOOLS and Mosek.
addpath(genpath('Path\to\SOSTOOLS'));
addpath(genpath('Path\to\MOSEK'));

% Number of Monte Carlo runs
N = 100;

% Create a timestamped subfolder to store results.
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
result_dir = fullfile(pwd, ['results_', timestamp]);
mkdir(result_dir);

% Initialize.
results_cell = cell(N, 1);
delta = zeros(N, 1);
valid = zeros(N, 1);

% Start parallel pool if not already running.
if isempty(gcp('nocreate'))
    parpool;
end

fprintf('### Started Monte Carlo simulation\n')
timer = tic;

% Run simulations.
parfor seed = 1:N
    try
        result = run_simulation(seed);
        results_cell{seed} = result;
        delta(seed) = result.delta;
        valid(seed) = result.valid;
    catch ME
        fprintf('Seed %d failed: %s\n', seed, ME.message);
        delta(seed) = NaN;
        valid(seed) = NaN;
    end
end

% Save results.
for seed = 1:N
    result = results_cell{seed};
    filename = fullfile(result_dir, sprintf('seed_%03d.mat', seed));
    save(filename, 'result');
end

% Get number of infeasible problems.
num_infeasible = sum(isnan(delta));
num_feasible = N-num_infeasible;

% Compute mean and variance of valid bounds.
lower_bound_perc = 100*(1-delta);
lower_bound_perc_mean = mean(lower_bound_perc, "omitnan");
lower_bound_perc_std = std(lower_bound_perc, "omitnan");
lower_bound_min = min(lower_bound_perc);
lower_bound_max = max(lower_bound_perc);

% Count number of valid barrier functions.
valid_barrier_certificates = sum(valid, "omitnan");

time_MC_simulation = toc(timer);
fprintf('### Monte Carlo simulation complete\nRuntime: %.2f s\n', time_MC_simulation);
fprintf('Number of times the SOS program is feasible: %i/%i (%.2f %%)\n', num_feasible, N, 100*(num_feasible/N))
fprintf('1-delta (mean ± std): %.2f %% ± %.2f %%\n', lower_bound_perc_mean, lower_bound_perc_std)
fprintf('1-delta range: [%.2f %%,  %.2f %%]\n', lower_bound_min, lower_bound_max)
fprintf('Number of times the certificate is valid for the true system: %i/%i (%.2f %%)\n', valid_barrier_certificates, N, 100*(valid_barrier_certificates/N))

% Save summary.
save(fullfile(result_dir, 'summary.mat'), 'N', 'delta', 'valid', 'num_infeasible', 'num_feasible', 'lower_bound_perc', 'lower_bound_perc_mean', 'lower_bound_perc_std', 'lower_bound_min', 'lower_bound_max', 'valid_barrier_certificates', 'time_MC_simulation');