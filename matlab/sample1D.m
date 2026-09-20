% SAMPLE1D - example of doing submovement decmoposition on 1D data
data = load('data/mirrorGame/tb_trial6.csv');

% data is sampled at 125 Hz
fs = 125;

% The first column in the x data - ignore the first 3 seconds
xdata = data(fs*3+1:end,1) ./ 59628 * 100; % convert to cm
time = (1:numel(xdata))'./fs;

% Apply a 4th order lowpass Butterworth filter with a cutoff of 4Hz
% This relatively low frequency is used to help remove tremor, etc
[B,A] = butter(2,4/(fs/2));
xdata_filtered = filtfilt(B,A,xdata);

% calculate the velocity (pad with a zero so it will keep the same length)
vel = [0;diff(xdata_filtered) ./ (1/fs)];

arng = [-60 60]; % cm/s
criteria = 0.03;
windowSize = 3; % seconds

% Optional fitting constraints / optimization settings
fittingConstraints = struct();

%% First fit the first 3 seconds
inds = time<=3;
submovementRange = 6:14;
[bestErrors,bestParameters,bestVelocity] = decompose1D(time(inds),vel(inds),submovementRange,arng,criteria);

%% Plot it
[~,ind] = min(bestErrors);
parameters = bestParameters{ind};
t0s = parameters(1:3:end-2);
Ds = parameters(2:3:end-1);
As = parameters(3:3:end);

figure;
h(1) = plot(time(inds),vel(inds));
hold on;
thist = (1:numel(bestVelocity{ind})) .* (time(2)-time(1));
h(2) = plot(thist,bestVelocity{ind})

for k=1:numel(t0s)
    thist = linspace(t0s(k),t0s(k)+Ds(k),100);
    thisv = minimumJerkVelocity1D(t0s(k),Ds(k),As(k),thist);
    h(3) = plot(thist,thisv,'Color',[0.5 0.5 0.5 0.5],'LineWidth',0.5);
end
legend('actual','fit','individual submovements');
%%
submovementRange = 6:14;
arng = [-60 60]; % cm/s
criteria = 0.03;
windowSize = 3; % seconds

[bestErrors,bestParameters,bestVelocity,decomposition] = decompose1Dwindows(time,vel,submovementRange,arng,criteria,windowSize,fittingConstraints);

%%
figure;
h(1) = plot(decomposition.time,decomposition.vel);
hold on;
h(2) = plot(decomposition.time,decomposition.reconstructedVelocity);
% Also plot the individual submovements
for k=1:numel(decomposition.t0s)
    ts = linspace(decomposition.t0s(k),decomposition.t0s(k)+decomposition.Ds(k),100);
    thisvel = minimumJerkVelocity1D(decomposition.t0s(k),...
        decomposition.Ds(k),...
        decomposition.As(k),ts);
    h(3) = plot(ts,thisvel,'Color',[0.5 0.5 0.5 0.5],'LineWidth',0.5);
end
xline(decomposition.startwindows,'Color',[0.5 0.5 0.5]);

legend(h,'Velocity','Reconstructed velocity','Submovements');
