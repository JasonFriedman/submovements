function [epsilon,grad,hess,sumpredicted,predictedx,predictedy] = calculateerrorMJ2D(parameters,time,vel,tangvel,timedelta)
% CALCULATEERRORMJ2D - calculate the error between the predicted and actual profile (in 2D)
% The predicted trajectory consists of the superposition of one or more minimum jerk velocity profiles
%
% [epsilon,grad,hess,sumpredicted,predictedx,predictedy] =  calculateerrorMJ2D(parameters,time,vel,tangvel,timedelta)
%
% The error is defined by (xvel - xpred)^2 + (yvel-ypred)^2 + (tangvel - tangpred)^2
%
% The function also optionally returns the gradient and Hessian
% (first-order and second-order partial derivatives), for use with
% optimization routines
%
% It can also optionally return the predicted minimum jerk trajectory
% (resulting from the superposition of the submovements)
%
% The parameters should be of length 4 * N (where N is the number submovements)
% each 4 parameters is T0 (onset time in seconds), D (duration in seconds),
% Ax (x amplitude) and Ay (y amplitude)
%
% time should be a 1 * N vector with the time of the recorded movement (in seconds)
%
% vel should be an N * 2 vector with the x and y velocities
%
% tangvel should contain the tangential velocity [i.e. tangvel = sqrt(vel(:,1).^2+vel(:,2).^2)  ]
%
% timedelta (optional, default = 0.005) is the time points to evaluate and
% compare the trajectories. It should match the time data [i.e. timedelta=time(2) - time(1)   ]

% Jason Friedman, 2021
% www.curiousjason.com

if nargin<5
    timedelta = 0.005; % 5 ms time delta
end

% Find end time of last submovement. If it is past end of the actual
% movement, pad actual movement with zeros

numsubmovements = length(parameters)/4;

lasttime = 0;
for k=1:numsubmovements
    % There are 4 parameters per submovement
    T0 = parameters(k*4-3);
    D =  parameters(k*4-2);
    lasttime = max([lasttime T0+D]);
end

% round lasttime to nearest 5ms
lasttime = round(lasttime* (1/timedelta))/ (1/timedelta);

if lasttime > time(end)
    time = [time(1:end-1); (time(end):timedelta:lasttime)'];
    vel(end+1:length(time),:) = 0;
    tangvel(end+1:length(time),:) = 0;
end

% Calculate the modelled trajectory

trajectoryx = vel(:,1);
trajectoryy = vel(:,2);

predictedx = zeros(numsubmovements,length(time));
predictedy = zeros(numsubmovements,length(time));
predicted = zeros(numsubmovements,length(time));

Jx= zeros(numsubmovements,4*numsubmovements,length(time));
Jy= zeros(numsubmovements,4*numsubmovements,length(time));
J= zeros(numsubmovements,4*numsubmovements,length(time));

Hx= zeros(numsubmovements,4*numsubmovements,4*numsubmovements,length(time));
Hy= zeros(numsubmovements,4*numsubmovements,4*numsubmovements,length(time));
H= zeros(numsubmovements,4*numsubmovements,4*numsubmovements,length(time));

for k=1:numsubmovements
    % There are 4 parameters per submovement
    T0 = parameters(k*4-3);
    D =  parameters(k*4-2);
    Dx = parameters(k*4-1);
    Dy = parameters(k*4);
    
    % find the appropriate time to calculate this over (T0 <= t <= T0+D)
    thisrng = find(time>T0 & time<T0+D);
    
    if nargout==1
        [predictedx(k,thisrng),predictedy(k,thisrng),predicted(k,thisrng)] ...
            = minimumJerkVelocity2D(T0,D,Dx,Dy,time(thisrng));
    elseif nargout==2
        [predictedx(k,thisrng),predictedy(k,thisrng),predicted(k,thisrng),...
            Jx(k,k*4-3:k*4,thisrng),Jy(k,k*4-3:k*4,thisrng),J(k,k*4-3:k*4,thisrng)] ...
            = minimumJerkVelocity2D(T0,D,Dx,Dy,time(thisrng));
    else
        [predictedx(k,thisrng),predictedy(k,thisrng),predicted(k,thisrng),...
            Jx(k,k*4-3:k*4,thisrng),Jy(k,k*4-3:k*4,thisrng),J(k,k*4-3:k*4,thisrng),...
            Hx(k,k*4-3:k*4,k*4-3:k*4,thisrng),Hy(k,k*4-3:k*4,k*4-3:k*4,thisrng),H(k,k*4-3:k*4,k*4-3:k*4,thisrng)] ...
            = minimumJerkVelocity2D(T0,D,Dx,Dy,time(thisrng));
    end
end
sumpredictedx = sum(predictedx,1)';
sumpredictedy = sum(predictedy,1)';
sumpredicted  = sum(predicted,1)';
sumtrajsq = sum(trajectoryx.^2 + trajectoryy.^2 + tangvel.^2);


if nargout>1
    sumJx = squeeze(sum(Jx,1));
    sumJy = squeeze(sum(Jy,1));
    sumJ = squeeze(sum(J,1));

    for k=1:size(sumJx,1)
        if size(sumpredictedx,1) ~= size(trajectoryx,1)
            keyboard;
        end
        grad(k,1) = 2/sumtrajsq * sum(...
            (sumpredictedx - trajectoryx).*sumJx(k,:)' + ...
            (sumpredictedy - trajectoryy).*sumJy(k,:)' + ...
            (sumpredicted - tangvel).*sumJ(k,:)');
    end
    if nargout>2
        sumHx = squeeze(sum(Hx,1));
        sumHy = squeeze(sum(Hy,1));
        sumH = squeeze(sum(H,1));
        for i=1:size(sumH,1)
            for j=1:size(sumH,2)
                hess(i,j) = 2/sumtrajsq * sum(...
                    sumJx(i,:).*sumJx(j,:) + ((sumpredictedx - trajectoryx).* squeeze(sumHx(i,j,:)))' + ...
                    sumJy(i,:).*sumJy(j,:) + ((sumpredictedy - trajectoryy).* squeeze(sumHy(i,j,:)))' + ...
                    sumJ(i,:) .*sumJ(j,:) +  ((sumpredicted - tangvel).* squeeze(sumH(i,j,:)))');
            end
        end
        
    end
end

epsilon = sum((sumpredictedx - trajectoryx).^2 + (sumpredictedy - trajectoryy).^2 ...
    +(sumpredicted - tangvel).^2) ./ sumtrajsq;

if nargout>3
    sumpredicted = [sumpredictedx sumpredictedy];
end