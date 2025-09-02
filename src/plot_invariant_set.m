function plot_invariant_set(h, f_true, r_min, r_max, varargin)
%PLOT_INVARIANT_SET Plot the forward invariant set X, the set X_min, the
% set X_max and a trajectory of the system.
%
%   Inputs:
%       h: barrier function
%       f_true: true state transition function
%       r_min: radius of the circular set X_min centered at the origin
%       r_max: radius of the circular set X_max centered at the origin
%
%   Variable-length input argument list:
%       x_trajectory: some state trajectory to be plotted
%       x_trajectory_label: label of the trajectory to be plotted

% Default values
x_trajectory = [];
x_trajectory_label = '';

% Read variable-length input argument list.
for i = 1:2:length(varargin)
    if strcmp('x_trajectory', varargin{i})
        x_trajectory = varargin{i+1};
    elseif strcmp('x_trajectory_label', varargin{i})
        x_trajectory_label = varargin{i+1};
    end
end

% Convert polynomial variables to matlab functions.
h = matlabFunction(pvar_to_sym(h));

% Define grid.
x = linspace(-r_max, r_max, 400);
y = linspace(-r_max, r_max, 400);
[X, Y] = meshgrid(x, y);

% Evaluate function on the grid.
Z = h(X, Y);

% Plot contour of h at level 0.
figure;
CM = contour(X, Y, Z, [0, 0], 'r', 'LineWidth', 1.5);
colormap([0.4660, 0.6740, 0.1880]);
hold on

% Plot X_min and X_max.
fill(r_min*sin(0:.1:2*pi), r_min*cos(0:.1:2*pi), [0.4660, 0.6740, 0.1880], FaceAlpha = .3)
fill(r_max*sin(0:.1:2*pi), r_max*cos(0:.1:2*pi), 'k', FaceAlpha = .3)

% Plot state trajectory if provided.
if ~isempty(x_trajectory)
    plot(x_trajectory(1, :), x_trajectory(2, :), 'Color', "#EDB120", LineWidth = 1.5)
end

% Compute system state one step into the future from selected points on 
% the boundary and add them to the plot.
pos = [];
i = 1;
while i < size(CM, 2)
    k = CM(2, i);
    for j = i + 1:i + k
        pos = [pos, CM(:, j)];
    end
    i = j + 1;
end
z = zeros(2, size(pos, 2));
zp = zeros(2, size(pos, 2));
for i = 1:5:size(pos, 2)
    x0 = pos(:, i);
    if (norm(x0) > r_max)
        continue
    end
    z(:, i) = x0;
    zp(:, i) = f_true(x0) - x0;
end

scale = 0.5;
zp = scale * zp;
quiver(z(1, :), z(2, :), zp(1, :), zp(2, :), 0, 'Color', "#D95319", HandleVisibility = 'off')

% Add legend etc.
legend('$\mathcal{C}$', '$\mathcal{X}_{min}$', '$\mathcal{X}_{max}$', x_trajectory_label, 'interpreter', 'latex');
xlabel('$x_1$', 'interpreter', 'latex')
ylabel('$x_2$', 'interpreter', 'latex')
