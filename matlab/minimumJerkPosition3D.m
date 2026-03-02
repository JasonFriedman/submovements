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
%    Az = displacement resulting from the movement (z)
%    x0 = starting position (x)
%    y0 = starting position (y)
%    z0 = starting position (z)
%
% The function is evaluated at times t
% Only value of t0 <= t <= t0 + D will be evaluated
% Earlier and later values will be set to x0 and x0 + Ax respectively for x
% and y0 and y0 + Ay respectively for y
% and z0 and z0 + Az respectively for z
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
[~,r,before,after,positionBasis] = minimumJerkBasis(t0,D,t);
x(r) = Ax * positionBasis(r) + x0;
y(r) = Ay * positionBasis(r) + y0;
z(r) = Az * positionBasis(r) + z0;

x(before) = x0; y(before) = y0; z(before) = z0;
x(after) = x0+Ax; y(after) = y0+Ay; z(after) = z0+Az;