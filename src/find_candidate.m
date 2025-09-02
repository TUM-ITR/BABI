function [h, q1, q2, q3, q4, q5, info] = find_candidate(MMH_samples_training, f_theta, n_x, r_min, r_max, gamma, orders)
%FIND_CANDIDATE Find barrier function candidate that is robust to all
% samples in the training set by solving an SOS program.
%
%   Inputs:
%       MMH_samples_training: samples in the training set
%       f_theta: state transition function parametrized by theta; has inputs (theta, x)
%       n_x: number of states
%       r_min: radius of the circular set X_min centered at the origin
%       r_max: radius of the circular set X_max centered at the origin
%       gamma: parameter in the class K function alpha(x) = gamma * x
%       orders: orders for the functions h, q_1, q_2, q_3, q_4, q_5
%
%   Outputs:
%       h: candidate barrier function
%       q1: function q_1
%       q2: function q_2
%       q3: function q_3
%       q4: function q_4
%       q5: function q_5
%       info: diagnostic information

%% Initialize
SOSP_timer = tic;
x = mpvar('x', n_x, 1);
prog = sosprogram(x);

%% Generate polynomials
[prog, h, coeff_h] = create_polynomial(prog, orders{1}, x);
[prog, q1, coeff_q1] = create_polynomial(prog, orders{2}, x);
[prog, q2, coeff_q2] = create_polynomial(prog, orders{2}, x);
[prog, q3, coeff_q3] = create_polynomial(prog, orders{2}, x);
[prog, q4, coeff_q4] = create_polynomial(prog, orders{3}, x(1));
[prog, q5, coeff_q5] = create_polynomial(prog, orders{3}, x(1));

%% Add constraints
relaxation1 = r_max^2 - (x' * x);
for k = 1:numel(MMH_samples_training)
    f = f_theta(MMH_samples_training{k}.theta, x);
    expr = subs(h, x, f) - (1 - gamma) * h;
    prog = sosineq(prog, expr-q1*relaxation1);
end

h_min = subs(q4, x(1), x'*x) - subs(q4, x(1), r_min^2);
h_max = subs(q5, x(1), x'*x) - subs(q5, x(1), r_max^2);
expr2 = h_max - h;
relaxation2 = r_min^2 - (x' * x);
expr3 = h - h_min - .001;

prog = sosineq(prog, expr2-q2*relaxation1);
prog = sosineq(prog, expr3-q3*relaxation2);
prog = sosineq(prog, q1);
prog = sosineq(prog, q2);
prog = sosineq(prog, q3);

% The following part constraints all parameters to lie between -10 and 10 for numerical stability.
for i = 1:length(coeff_h)
    prog = sosineq(prog, -coeff_h(i)+10);
    prog = sosineq(prog, coeff_h(i)+10);
end
for i = 1:length(coeff_q1)
    prog = sosineq(prog, -coeff_q1(i)+10);
    prog = sosineq(prog, coeff_q1(i)+10);
end
for i = 1:length(coeff_q2)
    prog = sosineq(prog, -coeff_q2(i)+10);
    prog = sosineq(prog, coeff_q2(i)+10);
end
for i = 1:length(coeff_q3)
    prog = sosineq(prog, -coeff_q3(i)+10);
    prog = sosineq(prog, coeff_q3(i)+10);
end
for i = 1:length(coeff_q4)
    prog = sosineq(prog, -coeff_q4(i)+10);
    prog = sosineq(prog, coeff_q4(i)+10);
end
for i = 1:length(coeff_q5)
    prog = sosineq(prog, -coeff_q5(i)+10);
    prog = sosineq(prog, coeff_q5(i)+10);
end

%% Solve SOS program
solver_opt.solver = 'mosek';
solver_opt.simplify = 'on';
solver_opt.params.MSK_IPAR_NUM_THREADS = 1;
solver_opt.params.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
solver_opt.params.MSK_DPAR_INTPNT_TOL_MU_RED = 1e-10;

fprintf('### Started SOS program\n')

prog = sossolve(prog, solver_opt);

time_SOSP = toc(SOSP_timer);
fprintf('### SOS program solved\nRuntime: %.2f s\n', time_SOSP);

%% Extract solution
h = sosgetsol(prog, h);
q1 = sosgetsol(prog, q1);
q2 = sosgetsol(prog, q2);
q3 = sosgetsol(prog, q3);
q4 = sosgetsol(prog, q4);
q5 = sosgetsol(prog, q5);
info = prog.solinfo.info;
end