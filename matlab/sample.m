%% Step 1: Download the source code and data from github and unzip it
websave('master.zip',...
    'https://github.com/JasonFriedman/submovements/archive/refs/heads/master.zip');
unzip('master.zip');
% Change to the directory
cd('submovements-master');
% Add the code directory to the path
addpath('matlab');

%% Step 2: load the data

[position,velocity,time] = loaddata('data/subject08day1pre/');

%% Step 3: Plot the data to become familiar with it
plotposition(position,time);
plotvelocity(velocity,time);

% Q: How many movements are there?
% Q: What are the units likely to be?
% Q: How long did the movements take (mean / std)

%% Step 4: Decompose a single trial into 4 submovements

numsubmovements = 4;

xrng = [-20 20];
yrng = [-10 10];

k=1;
[bestError1,bestParameters1,bestVelocity1] = decompose2D(time{k},velocity{k},numsubmovements,xrng,yrng);

% Q: how many parameters are there in total?
% Q: What are the units of the parameters?

%% Step 5: Plot the result

% velocity
figure;
plotSubmovements2D(bestParameters1,time{k},1);
hold on;
plot(time{k},velocity{k});

% position
figure;
x0 = position{k}(1,1);
y0 = position{k}(1,2);
plotSubmovements2D(bestParameters1,time{k},3,x0,y0);
hold on;
plot(time{k},position{k}(:,1),'b--','LineWidth',2);
plot(time{k},position{k}(:,2),'r--','LineWidth',2);

% position x vs y
figure;
plotSubmovements2D(bestParameters1,time{k},5,x0,y0);
hold on;
plot(position{k}(:,1),position{k}(:,2),'b--','LineWidth',2);

% Q: How well did the reconstruction match the original movement and does
% this differ for position and velocity. Why?

%% Step 6: Run the decompisition on all trials in the condition 
% (get a coffee!)
for k=1:numel(time)
    k
    [bestError(k),bestParameters(k,:),bestVelocity{k}] = decompose2D(time{k},velocity{k},numsubmovements,xrng,yrng);
end

%% Step 7: Plot the parameters as a function of time
parameterNames = {'t0','D','Ax','Ay'};
figure;
for submovement=1:4
    for parameter=1:4
        subplot(4,4,(submovement-1)*4+parameter);
        plot(bestParameters(:,(submovement-1)*4+parameter));
        if submovement==4
            xlabel('trial');
        end
        title([parameterNames{parameter} ' ' num2str(submovement)]);
    end
end

% Q: What trends can you observe in the data?

%% Step 8: Calculate two new parameters based on the decomposition: 
% overlap and relative onset time
%
% overlap: mean percent overlap between two submovements
overlapsPre = calculateOverlap(bestParameters); 

% relative onset time: time 2-4th submovement starts relative to duration
% of previous submovement

relativeOnsetsPre  = calculateRelativeOnset(bestParameters);

%% Step 9: Run the decomposition on the posttest and compare to the pretest
[positionPost,velocityPost,timePost] = loaddata('data/subject08day1post/');
% (get another coffee!)
for k=1:numel(timePost)
    k
    [bestErrorPost(k),bestParametersPost(k,:),bestVelocityPost{k}] = decompose2D(timePost{k},velocityPost{k},numsubmovements,xrng,yrng);
end

overlapsPost = calculateOverlap(bestParametersPost);
relativeOnsetsPost  = calculateRelativeOnset(bestParametersPost);

% Q: Is there a difference between the pre and post?


%% Step 10: Compare results to a subject from the observation group (34)
[positionPre2,velocityPre2,timePre2] = loaddata('data/subject34day1pre/');
% (get another coffee!)
for k=1:numel(timePre2)
    k
    [bestErrorPre2(k),bestParametersPre2(k,:),bestVelocityPre2{k}] = decompose2D(timePre2{k},velocityPre2{k},numsubmovements,xrng,yrng);
end

overlapsPre2 = calculateOverlap(bestParametersPre2);
relativeOnsetsPre2  = calculateRelativeOnset(bestParametersPre2);

%%

[positionPost2,velocityPost2,timePost2] = loaddata('data/subject34day1post/');
% (get another coffee!)
for k=1:numel(timePost2)
    k
    [bestErrorPost2(k),bestParametersPost2(k,:),bestVelocityPost2{k}] = decompose2D(timePost2{k},velocityPost2{k},numsubmovements,xrng,yrng);
end

overlapsPost2 = calculateOverlap(bestParametersPost2); 
relativeOnsetsPost2  = calculateRelativeOnset(bestParametersPost2);

% Q: For this subbject, is there a different between pre and post?
% What does this suggest about differences in coarticulation between
% the groups? 
% (Disclaimer: A proper analysis would of course look at all subjects
% and use an appropriate statistical test!)


