function [bestErrors,bestParameters,bestVelocity,decomposition] = decompose2Dwindows(time,vel,submovementRange,xrng,yrng,criteria,windowSize,fittingConstraints)
% DECOMPOSE2DWINDOWS - decompose two dimensional movement into submovements using velocity profiles
% divided into windows (useful for long duration movements)
%
% [bestErrors,bestParameters,bestVelocity,decomposition] = decompose2Dwindows(time,vel,submovementRange,xrng,yrng,criteria,windowSize)
%
% vel should be a N x 2 matrix, with the x and y velocities
%
% time should be a N x 1 vector with the corresponding time (in seconds)
%
% submovementRange is the number of submovements to look for, if it is
% empty or not specified, the function will try 1:4 submovements
%
% xrng is the valid range for the amplitude of x values (default = [-5 5])
%
% yrng is the valid range for the amplitude of y values (default = [0.1 5])
%
% min(t0) = 0.167 * submovement number
%
% criteria - stop if the bestError is less than the criteria (only relevant
% when numsubmovements is a vector of multiple values). This can save time
% by not checking for a higher number of submovements. Default is -inf
% (i.e. don't use), reasonable values could be 0.01-0.05
%
% windowSize - nominal duration of the window (in seconds) - default is 3 seconds
% The actual window progression is adaptive: if a fitted submovement extends
% beyond the current window, it is deferred to the next fit and the next
% window starts at that submovement onset
%
% bestErrors - the best (lowest) value of the error function (cell array - one
% per window)
%
% bestParameters contains the function parameters corresponding to the best values
% [t0 D Ax Ay]. If there are multiple submovements, it will have a
% length of 4*numsubmovements (cell array - one per window)
%
% bestVelocity is the velocity profile corresponding to the best values
% (cell array - one per window)
%
% decomposition - is a struct with a summary of the fits, with fields:
% t0s, Ds, Axs, Ays, endtimes, time, vel, startwindows, endwindows,
% submovementsVelocity, reconstructedVelocity

% Jason Friedman, 2026
% www.curiousjason.com

if nargin<3
    submovementRange = 1:4;
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

if nargin<7 || isempty(windowSize)
    windowSize = 3; % seconds
end

if nargin<8 || isempty(fittingConstraints)
    fittingConstraints = struct();
end

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2) ~= 2 || size(vel,1)==2
    error('velocity must be an N*2 matrix (it is a %d by %d matrix)',size(vel,1),size(vel,2));
end

if size(time,1) ~= size(vel,1)
    error('The length of the time vector (%d * %d) and the velocity matrix (%d * %d) must be equal',...
        size(time,1),size(time,2),size(vel,1),size(vel,2));
end

t0s = [];
Ds =  [];
Axs = [];
Ays = [];
endtimes = [];
startwindows = [];
endwindows = [];

tic
currentWindowStart = time(1);
w = 0;
while currentWindowStart < time(end)
    w = w+1;
    currentWindowEnd = min(currentWindowStart + windowSize,time(end));
    startwindows(w,1) = currentWindowStart;
    endwindows(w,1) = currentWindowEnd;

    thisinds = find(time>=currentWindowStart & time<=currentWindowEnd);
    if isempty(thisinds)
        break
    end

    thistime = time(thisinds);
    thisvel = vel(thisinds,:);

    % subtract any previously accepted submovements that overlap this window
    for k=1:numel(endtimes)
        if endtimes(k) > currentWindowStart
            [thisMJx,thisMJy] = minimumJerkVelocity2D(t0s(k),Ds(k),Axs(k),Ays(k),thistime);
            thisvel(:,1) = thisvel(:,1) - thisMJx';
            thisvel(:,2) = thisvel(:,2) - thisMJy';
        end
    end

    % decompose2D requires the time to start at zero, so subtract it here, then add it back on afterwards
    [windowErrors,windowParameters,windowVelocity] = decompose2D(thistime-thistime(1),thisvel,submovementRange,xrng,yrng,criteria,fittingConstraints);
    submovementInd = find(windowErrors<=criteria,1);
    if isempty(submovementInd)
        submovementInd = find(windowErrors<=0.05,1);
    end
    if isempty(submovementInd)
        submovementInd = find(windowErrors<0.1,1);
    end
    if isempty(submovementInd)
        [~,submovementInd] = min(windowErrors);
    end
    numSubmovements = submovementRange(submovementInd);

    % parameters are [t0 D Ax Ay]
    if ~isnan(numSubmovements) && ~isnan(windowParameters{submovementInd}(1))
        submovementParameters = reshape(windowParameters{submovementInd},4,numSubmovements)';
        thist0s = submovementParameters(:,1) + thistime(1); % put back in the right units
        thisDs = submovementParameters(:,2);
        thisAxs = submovementParameters(:,3);
        thisAys = submovementParameters(:,4);
        thisendtimes = thist0s+thisDs;

        keepMask = thisendtimes <= currentWindowEnd;
        deferredMask = ~keepMask;

        if any(deferredMask)
            nextWindowStart = min(thist0s(deferredMask));
            if nextWindowStart <= currentWindowStart
                % Avoid getting stuck when the only deferred submovement starts at the current boundary.
                keepMask(:) = true;
                deferredMask(:) = false;
                nextWindowStart = currentWindowEnd;
            else
                endwindows(w,1) = nextWindowStart;
            end
        else
            nextWindowStart = currentWindowEnd;
        end

        acceptedParameters = [thist0s(keepMask) thisDs(keepMask) thisAxs(keepMask) thisAys(keepMask)];
        acceptedLocalParameters = [acceptedParameters(:,1)-thistime(1) acceptedParameters(:,2) acceptedParameters(:,3) acceptedParameters(:,4)];
    else
        acceptedParameters = [];
        acceptedLocalParameters = [];
        nextWindowStart = currentWindowEnd;
    end

    if isempty(acceptedParameters)
        bestErrors{w} = NaN;
        bestParameters{w} = [];
        bestVelocity{w} = zeros(numel(thistime),2);
    else
        currentBestVelocity = zeros(numel(thistime),2);
        for k=1:size(acceptedLocalParameters,1)
            [thisVx,thisVy] = minimumJerkVelocity2D(...
                acceptedLocalParameters(k,1),acceptedLocalParameters(k,2),...
                acceptedLocalParameters(k,3),acceptedLocalParameters(k,4),...
                thistime-thistime(1));
            currentBestVelocity(:,1) = currentBestVelocity(:,1) + thisVx';
            currentBestVelocity(:,2) = currentBestVelocity(:,2) + thisVy';
        end
        currentBestError = sum(sum((currentBestVelocity - thisvel).^2,2)) / max(sum(sum(thisvel.^2,2)),1);

        bestErrors{w} = currentBestError;
        bestParameters{w} = reshape(acceptedParameters',1,[]);
        bestVelocity{w} = currentBestVelocity;

        t0s = [t0s;acceptedParameters(:,1)];
        Ds = [Ds;acceptedParameters(:,2)];
        Axs = [Axs;acceptedParameters(:,3)];
        Ays = [Ays;acceptedParameters(:,4)];
        endtimes = [endtimes;acceptedParameters(:,1)+acceptedParameters(:,2)];
    end

    endtime = toc;
    endtimehours = floor(endtime/60/60);
    endtime = endtime - endtimehours * 60 * 60;
    endtimeminutes = floor(endtime/60);
    endtime = endtime - endtimeminutes * 60;
    processedDuration = endwindows(w) - time(1);
    totalDuration = time(end) - time(1);
    processedPercent = processedDuration / totalDuration * 100;
    fprintf(['Finished window %d, time since start: %d hours, %d minutes, %d seconds, ' ...
        'processed %.1f%% (%.3f seconds from %.3f seconds)\n'],...
        w,endtimehours, endtimeminutes, round(endtime),...
        processedPercent,processedDuration,totalDuration);

    if nextWindowStart >= time(end) || processedDuration >= totalDuration-0.01
        break
    end
    currentWindowStart = nextWindowStart;
end

decomposition.t0s = t0s;
decomposition.Ds = Ds;
decomposition.Axs = Axs;
decomposition.Ays = Ays;
decomposition.endtimes = endtimes;
decomposition.time = time;
decomposition.vel = vel;
decomposition.startwindows = startwindows;
decomposition.endwindows = endwindows;

for k=numel(decomposition.t0s):-1:1
    [decomposition.submovementsVelocity(:,k,1),decomposition.submovementsVelocity(:,k,2)] = ...
        minimumJerkVelocity2D(...
        decomposition.t0s(k),decomposition.Ds(k),...
        decomposition.Axs(k),decomposition.Ays(k),...
        time);
end
decomposition.reconstructedVelocity = squeeze(sum(decomposition.submovementsVelocity,2));
