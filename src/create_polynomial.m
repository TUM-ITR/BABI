function [prog, polynomial, coefficients] = create_polynomial(prog, order, variables)
%WRITE_POLY Helper function that returns a polynomial of the specified order
% in the variables 'variables' along with its coefficients 'coefficients'.
%
%   Inputs:
%       prog: SOS program to which the polynomial is added
%       order: order of the polynomial
%       variables: variables that appear in the polynomial
%
%   Outputs:
%       prog: SOS program
%       polynomial: generated polynomial
%       coefficients: coefficients of the polynomial

Z = monomials(variables, order);
[prog, coefficients] = sospolymatrixvar(prog, monomials(variables, 0), size(Z));
polynomial = Z' * coefficients;
end
