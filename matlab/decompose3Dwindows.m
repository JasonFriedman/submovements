function [bestErrors,bestParameters,bestVelocity,decomposition] = decompose3Dwindows(time,vel,submovementRange,xrng,yrng,zrng,criteria,windowSize)
% DECOMPOSE3DWINDOWS - decompose three dimensional movement into submovements using the velocity profiles
% divided into windows (useful for long duration movements)
%
% [best,bestParameters,bestVelocity] = decompose3Dwindows(time,vel,submovementRange,xrng,yrng,zrng,criteria,windowSize)
%
% vel should be a N x 3 matrix, with the x, y and z velocities
%
% t should be a N x 1 matrix with the corresponding time (in seconds)
%
% submovementRange is the number of submovements to look for, if it is
% empty or not specified, the function will try 1 to 4 submovements
%
% xrng is the valid range for the amplitude of x values (default = [-5 5])
%
% yrng is the valid range for the amplitude of y values (default = [0.1 5])
%
% zrng is the valid range for the amplitude of y values (default = [-5 5])
%
% min(t0) = 0.167 * submovement number
%
% criteria - stop if the bestError is less than the criteria (only relevant
% when numsubmovements is a vector of multiple values). This can save time
% by not checking for a higher number of submovements
%
% window size - duration of the window (in seconds) - default is 3 seconds
%
% bestError the best (lowest) value of the error function (cell array - one
% per window)
%
% bestParameters contains the function parameters corresponding to the best values
% [t0 D Ax Ay Az]. If there are multiple submovements, it will be have a
% length of 5*numsubmovements (cell array - one per window)
%
% bestVelocity is the velocity profile coresponding to the best values
% (cell array - one per window)
%
% decomposition - is a struct with a summary of the fits, with fields:
% t0s, Ds, Axs, Ays, Azs, endtimes, time, vel, startwindows

% Jason Friedman, 2025
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

leftover = mod(time(end)-time(1),windowSize);

% if the last window is less than half the window size, just add it to the previous window
if leftover > 0.5
    numwindows = ceil( (time(end)-time(1))/windowSize);
    startwindows = 0:windowSize:windowSize*(numwindows-1);
    endwindows = startwindows + windowSize;
    endwindows(end) = time(end);
else
    numwindows = ceil( (time(end)-time(1))/windowSize)-1;
    startwindows = 0:windowSize:windowSize*(numwindows-1);
    endwindows = startwindows + windowSize;
    endwindows(end) = time(end);
end

t0s = [];
Ds =  [];
Axs = [];
Ays = [];
Azs = [];
endtimes = [];
thisendtimes = [];

tic
for w=1:numwindows
    thisinds = find(time>startwindows(w) & time<=endwindows(w));
    thistime = time(thisinds);
    thisvel = vel(thisinds,:);
    % subtract any submovements from the previous window
    for k=1:numel(thisendtimes)
        if thisendtimes(k) > startwindows(w)
            % subtract these from the current window
            [thisMJx,thisMJy,thisMJz] = minimumJerkVelocity3D(thist0s(k),thisDs(k),thisAxs(k),thisAys(k),thisAxs(k),thistime);
            thisvel(:,1) = thisvel(:,1) - thisMJx';
            thisvel(:,2) = thisvel(:,2) - thisMJy';
            thisvel(:,3) = thisvel(:,3) - thisMJz';
        end
    end

    % decompose3D requires the time to start at zero, so subtract it here, then add it back on afterwards
    [bestErrors{w},bestParameters{w},bestVelocity{w}] = decompose3D(thistime-thistime(1),thisvel,submovementRange,xrng,yrng,zrng,criteria);
    submovementInd = find(bestErrors{w}<=criteria,1);
    if isempty(submovementInd)
        submovementInd = find(bestErrors{w}<=0.05,1);
    end
    if isempty(submovementInd)
        submovementInd = find(bestErrors{w}<0.1,1);
    end
    if isempty(submovementInd)
        [~,submovementInd] = min(bestErrors{w});
    end
    numSubmovements = submovementRange(submovementInd);

    % parameters are [t0 D Ax Ay Az]
    submovementParameters = reshape(bestParameters{w}{submovementInd},5,numSubmovements)';
    thist0s = submovementParameters(:,1) + thistime(1); % put back in the right units
    thisDs = submovementParameters(:,2);
    thisAxs = submovementParameters(:,3);
    thisAys = submovementParameters(:,4);
    thisAzs = submovementParameters(:,5);
    thisendtimes = thist0s+thisDs;
    t0s      = [t0s;thist0s];
    Ds       = [Ds;thisDs];
    Axs      = [Axs;thisAxs];
    Ays      = [Ays;thisAys];
    Azs      = [Azs;thisAzs];
    endtimes = [endtimes;thisendtimes];
    endtime = toc;
    endtimehours = floor(endtime/60/60);
    endtime = endtime - endtimehours * 60 * 60;
    endtimeminutes = floor(endtime/60);
    endtime = endtime - endtimeminutes * 60;
    fprintf('Finished window %d of %d, time since start: %d hours, %d minutes, %d seconds\n',w,numwindows,endtimehours, endtimeminutes, round(endtime));
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