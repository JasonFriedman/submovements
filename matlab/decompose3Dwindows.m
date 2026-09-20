function [bestErrors,bestParameters,bestVelocity,decomposition] = decompose3Dwindows(time,vel,submovementRange,xrng,yrng,zrng,criteria,windowSize,fittingConstraints)
% DECOMPOSE3DWINDOWS - decompose three dimensional movement into submovements using the velocity profiles
% divided into windows (useful for long duration movements)
%
% [bestErrors,bestParameters,bestVelocity,decomposition] = decompose3Dwindows(time,vel,submovementRange,xrng,yrng,zrng,criteria,windowSize)
%
% vel should be a N x 3 matrix, with the x, y and z velocities
%
% time should be a N x 1 matrix with the corresponding time (in seconds)
%
% submovementRange is the number of submovements to look for, if it is
% empty or not specified, the function will try 1:4 submovements
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
% [t0 D Ax Ay Az]. If there are multiple submovements, it will have a
% length of 5*numsubmovements (cell array - one per window)
%
% bestVelocity is the velocity profile corresponding to the best values
% (cell array - one per window)
%
% decomposition - is a struct with a summary of the fits, with fields:
% t0s, Ds, Axs, Ays, Azs, endtimes, time, vel, startwindows, endwindows,
% submovementVelocity, reconstructedVelocity

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

if nargin<6 || isempty(zrng)
    zrng = [-5 5];
end

if nargin<7 || isempty(criteria)
    criteria = -inf;
end

if nargin<8 || isempty(windowSize)
    windowSize = 3; % seconds
end

if nargin<9 || isempty(fittingConstraints)
    fittingConstraints = struct();
end

if isfield(fittingConstraints,'minOnsetSpacing') && ~isempty(fittingConstraints.minOnsetSpacing)
    minOnsetSpacing = fittingConstraints.minOnsetSpacing;
else
    minOnsetSpacing = 0.167;
end
if minOnsetSpacing<=0
    error('fittingConstraints.minOnsetSpacing must be > 0');
end

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2) ~= 3 || size(vel,1)==3
    error('velocity must be an N*3 matrix (it is a %d by %d matrix)',size(vel,1),size(vel,2));
end

if size(time,1) ~= size(vel,1)
    error('The length of the time vector (%d * %d) and the velocity matrix (%d * %d) must be equal',...
        size(time,1),size(time,2),size(vel,1),size(vel,2));
end

t0s = [];
Ds =  [];
Axs = [];
Ays = [];
Azs = [];
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
            [thisMJx,thisMJy,thisMJz] = minimumJerkVelocity3D(t0s(k),Ds(k),Axs(k),Ays(k),Azs(k),thistime);
            thisvel(:,1) = thisvel(:,1) - thisMJx';
            thisvel(:,2) = thisvel(:,2) - thisMJy';
            thisvel(:,3) = thisvel(:,3) - thisMJz';
        end
    end

    % decompose3D requires the time to start at zero, so subtract it here, then add it back on afterwards
    [windowErrors,windowParameters,windowVelocity] = decompose3D(thistime-thistime(1),thisvel,submovementRange,xrng,yrng,zrng,criteria,fittingConstraints);
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

    % parameters are [t0 D Ax Ay Az]
    if ~isnan(numSubmovements) && ~isnan(windowParameters{submovementInd}(1))
        submovementParameters = reshape(windowParameters{submovementInd},5,numSubmovements)';
        thist0s = submovementParameters(:,1) + thistime(1); % put back in the right units
        thisDs = submovementParameters(:,2);
        thisAxs = submovementParameters(:,3);
        thisAys = submovementParameters(:,4);
        thisAzs = submovementParameters(:,5);
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

        acceptedParameters = [thist0s(keepMask) thisDs(keepMask) thisAxs(keepMask) thisAys(keepMask) thisAzs(keepMask)];

        % Enforce onset spacing against already accepted submovements from previous windows.
        if ~isempty(acceptedParameters)
            if isempty(t0s)
                previousOnset = -inf;
            else
                previousOnset = t0s(end);
            end
            spacingKeep = false(size(acceptedParameters,1),1);
            for kk=1:size(acceptedParameters,1)
                if acceptedParameters(kk,1) - previousOnset >= minOnsetSpacing-eps
                    spacingKeep(kk) = true;
                    previousOnset = acceptedParameters(kk,1);
                end
            end
            acceptedParameters = acceptedParameters(spacingKeep,:);
        end

        if isempty(acceptedParameters)
            acceptedLocalParameters = [];
        else
            acceptedLocalParameters = [acceptedParameters(:,1)-thistime(1) acceptedParameters(:,2) acceptedParameters(:,3) acceptedParameters(:,4) acceptedParameters(:,5)];
        end
    else
        acceptedParameters = [];
        acceptedLocalParameters = [];
        nextWindowStart = currentWindowEnd;
    end

    if isempty(acceptedParameters)
        bestErrors{w} = NaN;
        bestParameters{w} = [];
        bestVelocity{w} = zeros(numel(thistime),3);
    else
        currentBestVelocity = zeros(numel(thistime),3);
        for k=1:size(acceptedLocalParameters,1)
            [thisVx,thisVy,thisVz] = minimumJerkVelocity3D(...
                acceptedLocalParameters(k,1),acceptedLocalParameters(k,2),...
                acceptedLocalParameters(k,3),acceptedLocalParameters(k,4),acceptedLocalParameters(k,5),...
                thistime-thistime(1));
            currentBestVelocity(:,1) = currentBestVelocity(:,1) + thisVx';
            currentBestVelocity(:,2) = currentBestVelocity(:,2) + thisVy';
            currentBestVelocity(:,3) = currentBestVelocity(:,3) + thisVz';
        end
        currentBestError = sum(sum((currentBestVelocity - thisvel).^2,2)) / max(sum(sum(thisvel.^2,2)),1);

        bestErrors{w} = currentBestError;
        bestParameters{w} = reshape(acceptedParameters',1,[]);
        bestVelocity{w} = currentBestVelocity;

        t0s      = [t0s;acceptedParameters(:,1)];
        Ds       = [Ds;acceptedParameters(:,2)];
        Axs      = [Axs;acceptedParameters(:,3)];
        Ays      = [Ays;acceptedParameters(:,4)];
        Azs      = [Azs;acceptedParameters(:,5)];
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
decomposition.Azs = Azs;
decomposition.endtimes = endtimes;
decomposition.time = time;
decomposition.vel = vel;
decomposition.startwindows = startwindows;
decomposition.endwindows = endwindows;

for k=numel(decomposition.t0s):-1:1
    [decomposition.submovementsVelocity(:,k,1),decomposition.submovementsVelocity(:,k,2),...
        decomposition.submovementsVelocity(:,k,3)] = ...
        minimumJerkVelocity3D(...
        decomposition.t0s(k),decomposition.Ds(k),...
        decomposition.Axs(k),decomposition.Ays(k),...
        decomposition.Azs(k),...
        time);
end
decomposition.reconstructedVelocity = squeeze(sum(decomposition.submovementsVelocity,2));