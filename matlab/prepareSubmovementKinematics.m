function data = prepareSubmovementKinematics(parameters,dimensions,t,initialPositions)
% PREPARESUBMOVEMENTKINEMATICS - shared kinematic prep for submovement plotting

if dimensions ~= 2 && dimensions ~= 3
    error('dimensions must be 2 or 3');
end

parametersPerSubmovement = dimensions + 2;
if mod(numel(parameters),parametersPerSubmovement)~=0
    error('The parameters vector must have a length that is a multiple of %d',parametersPerSubmovement);
end

numSubmovements = numel(parameters)/parametersPerSubmovement;

parameterMatrix = reshape(parameters,parametersPerSubmovement,numSubmovements)';
t0 = parameterMatrix(:,1)';
D = parameterMatrix(:,2)';
A = parameterMatrix(:,3:end);

[~,order] = sort(t0);
t0 = t0(order);
D = D(order);
A = A(order,:);

starts = zeros(numSubmovements,dimensions);
starts(1,:) = initialPositions(:)';
for d=1:dimensions
    starts(2:numSubmovements,d) = starts(1,d) + cumsum(A(1:end-1,d));
end

tf = t0 + D;
if isempty(t)
    t = linspace(min(t0),max(tf),100);
end

vel = zeros(numSubmovements,numel(t),dimensions);
pos = zeros(numSubmovements,numel(t),dimensions);
for s = 1:numSubmovements
    if dimensions == 2
        [vx,vy] = minimumJerkVelocity2D(t0(s),D(s),A(s,1),A(s,2),t);
        [x,y] = minimumJerkPosition2D(t0(s),D(s),A(s,1),A(s,2),starts(s,1),starts(s,2),t);
        vel(s,:,1) = vx;
        vel(s,:,2) = vy;
        pos(s,:,1) = x;
        pos(s,:,2) = y;
    else
        [vx,vy,vz] = minimumJerkVelocity3D(t0(s),D(s),A(s,1),A(s,2),A(s,3),t);
        [x,y,z] = minimumJerkPosition3D(t0(s),D(s),A(s,1),A(s,2),A(s,3),starts(s,1),starts(s,2),starts(s,3),t);
        vel(s,:,1) = vx;
        vel(s,:,2) = vy;
        vel(s,:,3) = vz;
        pos(s,:,1) = x;
        pos(s,:,2) = y;
        pos(s,:,3) = z;
    end
end

posRelative = pos;
for s=2:numSubmovements
    for d=1:dimensions
        posRelative(s,:,d) = pos(s,:,d) - starts(s,d);
    end
end

data.numSubmovements = numSubmovements;
data.t0 = t0;
data.D = D;
data.A = A;
data.starts = starts;
data.t = t;
data.vel = vel;
data.pos = pos;
data.sumPos = squeeze(sum(posRelative,1));
if dimensions == 2
    data.sumPos = reshape(data.sumPos, numel(t), 2);
else
    data.sumPos = reshape(data.sumPos, numel(t), 3);
end
