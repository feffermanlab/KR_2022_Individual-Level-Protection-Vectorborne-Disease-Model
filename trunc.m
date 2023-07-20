function y = trunc(x, n)
% Truncate matrix/scalar x to n decimal places
if nargin < 2; n = 0; end % default to standard fix behaviour if no n given
y = fix(x*10^n)/10^n;      % return value truncated to n decimal places
end