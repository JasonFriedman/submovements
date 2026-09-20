function [epsilon,grad,hess,sumpredicted,predicted] = calculateerrorMJ1D(parameters,time,vel,timedelta)
% CALCULATEERRORMJ1D - calculate the error between the predicted and actual profile (in 1D)
% The predicted trajectory consists of the superposition of one or more minimum jerk velocity profiles
%
% [epsilon,grad,hess,sumpredicted,predicted] =  calculateerrorMJ1D(parameters,time,vel,timedelta)
%
% The error is defined by (vel - pred)^2 + (abs(vel) - abs(pred))^2
%
% The function also optionally returns the gradient and Hessian
% (first-order and second-order partial derivatives), for use with
% optimization routines
%
% It can also optionally return the predicted minimum jerk trajectory
% (resulting from the superposition of the submovements)
%
% The parameters should be of length 3 * N (where N is the number submovements)
% each 3 parameters is T0 (onset time in seconds), D (duration in seconds),
% and A (displacement)
%
% time should be a N * 1 vector with the time of the recorded movement (in seconds)
%
% vel should be an N * 1 vector with the velocity
%
% timedelta (optional, default = 0.005) is the time points to evaluate and
% compare the trajectories. It should match the time data [i.e. timedelta=time(2) - time(1)]
%
% Jason Friedman, 2026
% www.curiousjason.com

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
T = length(time);
np = 3*numsubmovements;

predicted = zeros(numsubmovements,T);

% sumJ/sumH are built directly in their final (block-diagonal) shape -
% each submovement only ever contributes to its own 3 columns/rows, so
% the extra "numsubmovements" leading dimension used previously (and then
% summed away) just wasted an N-fold amount of memory and computation.
if nargout>1
    sumJ = zeros(np,T);
end
if nargout>2
    sumH = zeros(np,np,T);
    sumabsH = zeros(np,np,T);
end

for k=1:numsubmovements
    % There are 3 parameters per submovement
    T0 = parameters(k*3-2);
    D =  parameters(k*3-1);
    A =  parameters(k*3);
    cols = k*3-2:k*3;

    % find the appropriate time to calculate this over (T0 <= t <= T0+D)
    thisrng = find(time>=T0 & time<=T0+D);

    if nargout==1
        predicted(k,thisrng) = minimumJerkVelocity1D(T0,D,A,time(thisrng));
    elseif nargout==2
        [predicted(k,thisrng),Jk] = minimumJerkVelocity1D(T0,D,A,time(thisrng));
        sumJ(cols,thisrng) = Jk;
    else
        [predicted(k,thisrng),Jk,Hk] = minimumJerkVelocity1D(T0,D,A,time(thisrng));
        sumJ(cols,thisrng) = Jk;
        sumH(cols,cols,thisrng) = Hk;
        signk = sign(predicted(k,thisrng));
        sumabsH(cols,cols,thisrng) = Hk .* reshape(signk,[1 1 numel(thisrng)]);
    end
end

sumpredicted = sum(predicted,1)';
sumabspredicted = sum(abs(predicted),1)';
sumtrajsq  = sum(trajectory.^2);
if sumtrajsq==0
    sumtrajsq = 1;
end

if nargout>1
    signpredicted = sign(predicted);
    sumabsJ = zeros(np,T);
    for k=1:numsubmovements
        cols = k*3-2:k*3;
        sumabsJ(cols,:) = sumJ(cols,:) .* signpredicted(k,:);
    end

    errTerm = sumpredicted - trajectory;
    absErrTerm = sumabspredicted - abs(trajectory);

    grad = zeros(np,1);
    for k=1:np
        grad(k,1) = 2/sumtrajsq * sum(...
            errTerm.*sumJ(k,:)' + ...
            absErrTerm.*sumabsJ(k,:)');
    end

    if nargout>2
        hess = zeros(np,np);
        for i=1:np
            for j=1:np
                hess(i,j) = 2/sumtrajsq * sum(...
                    sumJ(i,:).*sumJ(j,:) + (errTerm.* squeeze(sumH(i,j,:)))' + ...
                    sumabsJ(i,:).*sumabsJ(j,:) + (absErrTerm.* squeeze(sumabsH(i,j,:)))');
            end
        end
    end
end

epsilon = sum((sumpredicted - trajectory).^2 + ...
    (sumabspredicted-abs(trajectory)).^2) ./ sumtrajsq;
end
