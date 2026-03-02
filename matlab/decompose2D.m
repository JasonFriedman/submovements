function [bestError,bestParameters,bestVelocity] = decompose2D(time,vel,numsubmovements,xrng,yrng,criteria)
% DECOMPOSE2D - decompose two dimensional movement into submovements using the velocity profiles
%
% [bestError,bestParameters,bestVelocity] = decompose2D(time,vel,numsubmovements,xrng,yrng,criteria)
%
% vel should be a N x 2 matrix, with the x and y velocities
%
% time should be a N x 1 matrix with the corresponding time (in seconds)
%
% numsubmovements is the number of submovements to look for, if it is
% empty or not specified, the function will try 1 to 4 submovements
%
% xrng is the valid range for the amplitude of x values (default = [-5 5])
%
% yrng is the valid range for the amplitude of y values (default = [0.1 5])
%
% min(t0) = 0.167 * submovement number
%
% criteria - stop if the bestError is less than the criteria (only relevant
% when numsubmovements is a vector of multiple values). This can save time
% by not checking for a higher number of submovements
%
%
% bestError the best (lowest) value of the error function
%
% bestParameters contains the function parameters corresponding to the best values
% [t0 D Ax Ay]. If there are multiple submovements, it will have a
% length of 4*numsubmovements
%
% bestVelocity is the velocity profile corresponding to the best values

% Jason Friedman, 2021
% www.curiousjason.com

if nargin<3
    numsubmovements = [];
end

if nargin<4 || isempty(xrng)
    xrng = [-5 5];
end

if nargin<5 || isempty(yrng)
    yrng = [0.1 5];
end

if nargin<6 || isempty(criteria)
    criteria = -inf;
end

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2) >2 || size(vel,1)==2
    error('vel must be an N*2 matrix');
end

if size(time,1) ~= size(vel,1)
    error('The length of the time vector and the vel matrix must be equal');
end 

if isempty(numsubmovements) || length(numsubmovements)>1
    % If not specified, try 1 to 4 submovements
    if isempty(numsubmovements)
        numsubmovements = 1:4;
    end
    bestError = NaN * ones(1,max(numsubmovements)); bestParameters = cell(1,max(numsubmovements)); bestVelocity = cell(1,max(numsubmovements));
    
    for k=1:numel(numsubmovements)
        [bestError(k),bestParameters{k},bestVelocity{k}] = decompose2D(time,vel,numsubmovements(k),xrng,yrng);
        if bestError(k)<criteria
            return
        end
    end
    return
end

% parameters are T0, D, Ax Ay
% ranges are
% 0 <= T0 <= finaltime - 0.167
% 0.167 <= D <= finaltime
% xrng(1) <= Ax <= xrng(2)
% yrng(1) <= Ay <= yrng(2)
lb_0 = [0                           0.167     xrng(1) yrng(1)];
ub_0 = [max([time(end)-0.167 0.1])  1.0       xrng(2) yrng(2)];
pps = 4; % parameters per submovement

v = vel(:,1:2);
tv = sqrt(vel(:,1).^2 + vel(:,2).^2);

[bestError,bestParameters,bestVelocity] = decomposeND(...
    time,numsubmovements,lb_0,ub_0,pps,...
    @(parameters) calculateerrorMJ2D(parameters,time,v,tv));