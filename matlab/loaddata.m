function [positionfiltered,velocity,time] = loaddata(dirname)

d = dir([dirname '/*.csv']);

if isempty(d) 
    error('Must specify a directory to load the csv files from');
end

clear block trial
for k=1:numel(d)
   filename{k} = d(k).name; 
   r = regexp(filename{k},'tb_.*block(\d*)_trial(\d*).csv','tokens');
   block(k) = str2double(r{1}{1}); 
   trial(k) = str2double(r{1}{2});
end

for b=1:max(block)
    for t=1:max(trial)
        trialindex = find(b==block & t==trial);
        if isempty(trialindex)
            continue;
        end
        trialnum = (b-1)*max(trial) + t;
        data = load([d(trialindex).folder '/' d(trialindex).name]);
        pressure = data(:,4);
        % Only take the part where the pressure > 0
        position{trialnum} = data(pressure>0,1:2)./1000;
        % Calculate velocity after filtering
        
        time{trialnum} = data(pressure>0,5) ./ 1000; % seconds
        time{trialnum} = time{trialnum} - time{trialnum}(1);
        dt = median(diff(time{trialnum}));
        [B,A] = butter(2,5/((1/dt)/2));
        positionfiltered{trialnum} = filtfilt(B,A,position{trialnum});
        velocity{trialnum} = [[0 0]; diff(positionfiltered{trialnum})./dt];
    end
end
    