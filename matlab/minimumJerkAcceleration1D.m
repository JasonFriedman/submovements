function a = minimumJerkAcceleration1D(t0,D,A,t)
% minimumJerkAcceleration1D - evaluate the acceleration for a 1D minimum jerk curve
%
% a = minimumJerkAcceleration1D(t0,D,A,t)
%
% see Flash and Hogan (1985) for details on the minimum jerk equation
%
%    t0 = movement start time
%    D  = movement duration
%    Ax = displacement resulting from the movement
%
% The function is evaluated at times t
% Only value of t0 < t < t0 + D will be evaluated (rest will be zero)

% Jason Friedman, 2024
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


% nt is normalized time (0 <= nt <= 1)
a = zeros(1,numel(t));
nt = (t-t0)./D;
r = (nt>=0 & nt<=1);
a(r) = A/(D^2) * (-180 * nt(r).^2 + 120 * nt(r).^3 + 60 * nt(r));