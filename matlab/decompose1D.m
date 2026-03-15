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

function [epsilon,grad,hess,sumpredicted,predicted] = calculateerrorMJ1D(parameters,time,vel,timedelta)
% CALCULATEERRORMJ1D - calculate 1D velocity reconstruction error for minimum jerk submovements

if nargin<4
    timedelta = 0.005; % 5 ms time delta
end

numsubmovements = length(parameters)/3;

lasttime = 0;
for k=1:numsubmovements
    % There are 3 parameters per submovement
    T0 = parameters(k*3-2);
    D =  parameters(k*3-1);
    lasttime = max([lasttime T0+D]);
end

% round lasttime to nearest timedelta
lasttime = round(lasttime* (1/timedelta))/ (1/timedelta);

if lasttime > time(end)
    time = [time(1:end-1); (time(end):timedelta:lasttime)'];
    vel(end+1:length(time),:) = 0;
end

trajectory = vel(:,1);

predicted = zeros(numsubmovements,length(time));

J = zeros(numsubmovements,3*numsubmovements,length(time));
H = zeros(numsubmovements,3*numsubmovements,3*numsubmovements,length(time));

for k=1:numsubmovements
    % There are 3 parameters per submovement
    T0 = parameters(k*3-2);
    D =  parameters(k*3-1);
    A =  parameters(k*3);

    % find the appropriate time to calculate this over (T0 <= t <= T0+D)
    thisrng = find(time>=T0 & time<=T0+D);

    if nargout==1
        predicted(k,thisrng) = minimumJerkVelocity1D(T0,D,A,time(thisrng));
    elseif nargout==2
        [predicted(k,thisrng),J(k,k*3-2:k*3,thisrng)] = minimumJerkVelocity1D(T0,D,A,time(thisrng));
    else
        [predicted(k,thisrng),J(k,k*3-2:k*3,thisrng),H(k,k*3-2:k*3,k*3-2:k*3,thisrng)] = minimumJerkVelocity1D(T0,D,A,time(thisrng));
    end
end

sumpredicted = sum(predicted,1)';
sumtrajsq = sum(trajectory.^2);
if sumtrajsq==0
    sumtrajsq = 1;
end

if nargout>1
    sumJ = squeeze(sum(J,1));

    for k=1:size(sumJ,1)
        grad(k,1) = 2/sumtrajsq * sum((sumpredicted - trajectory).*sumJ(k,:)');
    end

    if nargout>2
        sumH = squeeze(sum(H,1));
        for i=1:size(sumH,1)
            for j=1:size(sumH,2)
                hess(i,j) = 2/sumtrajsq * sum(...
                    sumJ(i,:).*sumJ(j,:) + ((sumpredicted - trajectory).* squeeze(sumH(i,j,:)))');
            end
        end
    end
end

epsilon = sum((sumpredicted - trajectory).^2) ./ sumtrajsq;

end
