function symfun = pvar_to_sym(dpvarfun)
%PVAR_TO_SYM Helper function that converts a polynomial to a symbolic
% expression.
%
%   Inputs:
%       dpvarfun: polynomial (struct)
%
%   Outputs:
%       symfun: symbolic expression

coeff = dpvarfun.coefficient;
deg = dpvarfun.degmat;
vars = dpvarfun.varname;
x = sym("x", [numel(vars), 1], 'real');
monomial_vector = 1;

for i = 1:numel(vars)
    monomial_vector = monomial_vector .* x(i).^deg(:, i);
end

symfun = monomial_vector' * coeff;
end
