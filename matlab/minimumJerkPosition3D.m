function [x,y,z] = minimumJerkPosition3D(t0,D,Ax,Ay,Az,x0,y0,z0,t)
% minimumJerkPosition3D - evaluate a 3D minimum jerk curve
%
% [x,y,z] = minimumJerkPosition3D(t0,D,Ax,Ay,Az,x0,y0,z0,t)
%
% see Flash and Hogan (1985) for details on the minimum jerk equation
%
%    t0 = movement start time
%    D  = movement duration
%    Ax = displacement resulting from the movement (x)
%    Ay = displacement resulting from the movement (y)
%    Az = displacement resulting from the movement (y)
%    x0 = starting position (x)
%    y0 = starting position (y)
%    z0 = starting position (y)
%
% The function is evaluated at times t
% Only value of t0 < t < t0 + D will be evaluated (rest will be zero)
%

% Jason Friedman, 2021
% www.curiousjason.com

if numel(t0)>1
    error('t0 must be a number, not a matrix');
end

if numel(D)>1
    error('D must be a number, not a matrix');
end

if numel(Ax)>1
    error('Ax must be a number, not a matrix');
end

if numel(Ay)>1
    error('Ay must be a number, not a matrix');
end

if numel(Az)>1
    error('Az must be a number, not a matrix');
end

if numel(x0)>1
    error('x0 must be a number, not a matrix');
end

if numel(y0)>1
    error('y0 must be a number, not a matrix');
end

if numel(z0)>1
    error('z0 must be a number, not a matrix');
end

% nt is normalized time (0 <= nt <= 1)
x = zeros(1,numel(t));
nt = (t-t0)./D;
r = (nt>=0 & nt<=1);
before = nt<0;
after = nt>1;
x(r) = Ax * (-15 * nt(r).^4 + 6 * nt(r).^5 + 10 * nt(r).^3) + x0;
y(r) = Ay * (-15 * nt(r).^4 + 6 * nt(r).^5 + 10 * nt(r).^3) + y0;
z(r) = Az * (-15 * nt(r).^4 + 6 * nt(r).^5 + 10 * nt(r).^3) + z0;

x(before) = x0; y(before) = y0; z(before) = z0;
x(after) = x0+Ax; y(after) = y0+Ay; z(after) = z0+Az;