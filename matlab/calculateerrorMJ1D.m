function [epsilon,grad,hess,sumpredicted,predicted] = calculateerrorMJ1D(parameters,time,vel,timedelta)
% CALCULATEERRORMJ1D - calculate the error between the predicted and actual profile (in 1D)
% The predicted trajectory consists of the superposition of one or more minimum jerk velocity profiles
%
% [epsilon,grad,hess,sumpredicted,predicted] =  calculateerrorMJ1D(parameters,time,vel,timedelta)
%
% The error is defined by (vel - pred)^2
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
