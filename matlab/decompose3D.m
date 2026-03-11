function [bestError,bestParameters,bestVelocity] = decompose3D(time,vel,numsubmovements,xrng,yrng,zrng,criteria,fittingConstraints)
% DECOMPOSE3D - decompose three dimensional movement into submovements using the velocity profiles
%
% [bestError,bestParameters,bestVelocity] = decompose3D(time,vel,numsubmovements,xrng,yrng,zrng,criteria)
%
% vel should be a N x 3 matrix, with the x, y and z velocities
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
% zrng is the valid range for the amplitude of z values (default = [-5 5])
%
% min(t0) = 0.167 * submovement number
%
% criteria - stop if the bestError is less than the criteria (only relevant
% when numsubmovements is a vector of multiple values). This can save time
% by not checking for a higher number of submovements
%
% bestError the best (lowest) value of the error function
%
% bestParameters contains the function parameters corresponding to the best values
% [t0 D Ax Ay Az]. If there are multiple submovements, it will have a
% length of 5*numsubmovements
%
% bestVelocity is the velocity profile corresponding to the best values

% Jason Friedman, 2024
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

if nargin<6 || isempty(zrng)
    zrng = [-5 5];
end

if nargin<7 || isempty(criteria)
    criteria = -inf;
end

if nargin<8 || isempty(fittingConstraints)
    fittingConstraints = struct();
end

constraints = resolveFittingConstraints(fittingConstraints);

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2) > 3 || size(vel,1)==3
    error('vel must be an N*3 matrix');
end

if size(time,1) ~= size(vel,1)
    error('The length of the time vector and the vel matrix must be equal');
end 

if isempty(numsubmovements) || length(numsubmovements)>1
    % If not specified, try 1 to 4 submovements
    if isempty(numsubmovements)
        numsubmovements = 1:4;
    end
    bestError = NaN * ones(1,numel(numsubmovements)); bestParameters = cell(1,numel(numsubmovements)); bestVelocity = cell(1,numel(numsubmovements));
    
    for k=1:numel(numsubmovements)
        [bestError(k),bestParameters{k},bestVelocity{k}] = decompose3D(time,vel,numsubmovements(k),xrng,yrng,zrng,criteria,fittingConstraints);
        if bestError(k)<criteria
            return
        end
    end
    return
end

% parameters are T0, D, Ax Ay, Az
% ranges are
% 0 <= T0 <= finaltime - minOnsetSpacing
% minDuration <= D <= finaltime
% xrng(1) <= Ax <= xrng(2)
% yrng(1) <= Ay <= yrng(2)
% zrng(1) <= Az <= zrng(2)
lb_0 = [0                                                            constraints.minDuration      xrng(1) yrng(1) zrng(1)];
ub_0 = [max([time(end)-constraints.minOnsetSpacing constraints.minUpperBoundTime])  constraints.maxDuration  xrng(2) yrng(2) zrng(2)];
pps = 5; % parameters per submovement

if numel(time)==0
    bestError = NaN;
    bestParameters = NaN(1,numsubmovements*5);
    bestVelocity = NaN;
    return
end

v = vel(:,1:3);
tv = sqrt(vel(:,1).^2 + vel(:,2).^2 + vel(:,3).^2);
timedelta = time(2)-time(1);

[bestError,bestParameters,bestVelocity] = decomposeND(...
    time,numsubmovements,lb_0,ub_0,pps,...
    @(parameters) calculateerrorMJ3D(parameters,time,v,tv,timedelta),fittingConstraints);