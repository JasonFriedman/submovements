function [bestError,bestParameters,bestVelocity] = decompose1D(time,vel,numsubmovements,arng,criteria,fittingConstraints)
% DECOMPOSE1D - decompose one dimensional movement into submovements using the velocity profile
%
% [bestError,bestParameters,bestVelocity] = decompose1D(time,vel,numsubmovements,arng,criteria)
%
% vel should be a N x 1 vector, with the velocity profile
%
% time should be a N x 1 vector with the corresponding time (in seconds)
%
% numsubmovements is the number of submovements to look for, if it is
% empty or not specified, the function will try 1 to 4 submovements
%
% arng is the valid range for the amplitude values (default = [-5 5])
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
% [t0 D A]. If there are multiple submovements, it will have a
% length of 3*numsubmovements
%
% bestVelocity is the velocity profile corresponding to the best values

% Jason Friedman, 2026
% www.curiousjason.com

if nargin<3
    numsubmovements = [];
end

if nargin<4 || isempty(arng)
    arng = [-5 5];
end

if nargin<5 || isempty(criteria)
    criteria = -inf;
end

if nargin<6 || isempty(fittingConstraints)
    fittingConstraints = struct();
end

constraints = resolveFittingConstraints(fittingConstraints);

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2)>1
    error('vel must be an N*1 vector');
end

if size(time,1) ~= size(vel,1)
    error('The length of the time vector and the vel vector must be equal');
end

if isempty(numsubmovements) || length(numsubmovements)>1
    % If not specified, try 1 to 4 submovements
    if isempty(numsubmovements)
        numsubmovements = 1:4;
    end
    bestError = NaN * ones(1,numel(numsubmovements));
    bestParameters = cell(1,numel(numsubmovements));
    bestVelocity = cell(1,numel(numsubmovements));

    for k=1:numel(numsubmovements)
        [bestError(k),bestParameters{k},bestVelocity{k}] = decompose1D(time,vel,numsubmovements(k),arng,criteria,fittingConstraints);
        if bestError(k)<criteria
            return
        end
    end
    return
end

% parameters are T0, D, A
% ranges are
% 0 <= T0 <= finaltime - minOnsetSpacing
% minDuration <= D <= finaltime
% arng(1) <= A <= arng(2)
lb_0 = [0                                                            constraints.minDuration      arng(1)];
ub_0 = [max([time(end)-constraints.minOnsetSpacing constraints.minUpperBoundTime])  constraints.maxDuration  arng(2)];
pps = 3; % parameters per submovement

v = vel(:,1);
timedelta = time(2)-time(1);

[bestError,bestParameters,bestVelocity] = decomposeND(...
    time,numsubmovements,lb_0,ub_0,pps,...
    @(parameters) calculateerrorMJ1D(parameters,time,v,timedelta),fittingConstraints);

end
