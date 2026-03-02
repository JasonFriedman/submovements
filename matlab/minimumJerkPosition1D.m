function x = minimumJerkPosition1D(t0,D,A,x0,t)
% minimumJerkPosition1D - evaluate a 1D minimum jerk curve
%
% x = minimumJerkPosition1D(t0,D,A,x0,t)
%
% see Flash and Hogan (1985) for details on the minimum jerk equation
%
%    t0 = movement start time
%    D  = movement duration
%    A = displacement resulting from the movement
%    x0 = starting position
%
% The function is evaluated at times t
% Only value of t0 <= t <= t0 + D will be evaluated 
% Earlier and later values will be set to x0 and x0 + A respectively
%

% Jason Friedman, 2021
% www.curiousjason.com

if numel(t0)>1
    error('t0 must be a number, not a matrix');
end

if numel(D)>1
    error('D must be a number, not a matrix');
end

if numel(A)>1
    error('A must be a number, not a matrix');
end

if numel(x0)>1
    error('x0 must be a number, not a matrix');
end

% nt is normalized time (0 <= nt <= 1)
x = zeros(1,numel(t));
[~,r,before,after,positionBasis] = minimumJerkBasis(t0,D,t);
x(r) = A * positionBasis(r) + x0;
x(before) = x0;
x(after) = x0+A;

