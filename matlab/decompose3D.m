function [bestError,bestParameters,bestVelocity] = decompose3D(time,vel,numsubmovements,xrng,yrng,zrng)
% DECOMPOSE - decompose three dimensional movement into submovements using the velocity profiles
%
% [best,bestParameters,bestVelocity] = decompose(time,vel,numsubmovements,xrng,yrng,zrng)
%
% vel should be a N x 3 matrix, with the x, y and z velocities
%
% t should be a N x 1 matrix with the corresponding time (in seconds)
%
% numsubmovements is the number of submovements to look for, if it is
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
%
% bestError the best (lowest) value of the error function
%
% bestParameters contains the function parameters corresponding to the best values
% [t0 D Ax Ay Az]. If there are multiple submovements, it will be have a
% length of 5*numsubmovements
%
% bestVelocity is the velocity profile coresponding to the best values

% Jason Friedman, 2023
% www.curiousjason.com

if nargin<3
    numsubmovements = [];
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

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2) > 3 || size(vel,1)==3
    error('vel must be an N*3 matrix');
end

if size(time,1) ~= size(vel,1)
    error('The length of the time vector and the vel matrix must be equal');
end 

if isempty(numsubmovements) || length(numsubmovements)>1
    % If not specified, try 1 to 4 submovements
    if isempty(numsubmovements)
        numsubmovements = 1:4;
    end
    bestError = NaN * ones(1,max(numsubmovements)); bestParameters = cell(1,max(numsubmovements)); bestVelocity = cell(1,max(numsubmovements));
    
    for k=numsubmovements
        [bestError(k),bestParameters{k},bestVelocity{k}] = decompose3D(time,vel,k,xrng,yrng,zrng);
    end
    return;
end

if numel(time)==0
    bestError = NaN;
    bestParameters = NaN(1,numsubmovements*5);
    bestVelocity = NaN;
    return
end
    

bestError = inf;
ignoreerrors = true;

% parameters are T0, D, Ax Ay, Az
% ranges are
% 0 <= T0 <= finaltime - 0.167
% 0.167 <= D <= finaltime
% xrng(1) <= Ax <= xrng(2)
% yrng(1) <= Ay <= yrng(2)
% zrng(1) <= Az <= zrng(2)
lb_0 = [0                           0.167     xrng(1) yrng(1) zrng(1)];
ub_0 = [max([time(end)-0.167 0.1])  1.0       xrng(2) yrng(2) zrng(2)];
pps = 5; % parameters per submovement

if any(lb_0>ub_0)
    error('Lower bounds exceed upper bound - infeasible');
end

% In Roher & Hogan 2006, they selected 10 random parameter selections
% Here we use 20 (increases a lot the likelihood to converge to the same
% solution on multiple runs)
count=1;
while count<=20
    for i=1:numsubmovements
        % Randomly select 20 starting positions in the legal range
        initialparameters(1,i*pps-(pps-1):i*pps) = lb_0 + (ub_0-lb_0) .* rand(1,pps);
        thislb_0 = lb_0;
        thislb_0(1) = (i-1)*0.167;
        lb(i*pps-(pps-1):i*pps) = thislb_0;
        ub(i*pps-(pps-1):i*pps) = ub_0;
    end
    v = vel(:,1:3);
    tv = sqrt(vel(:,1).^2 + vel(:,2).^2 + vel(:,3).^2);
    % Turn on GradObj and Hessian to use the gradient and Hessian
    % calculated in the functions
    options = optimset('GradObj','on','Hessian','on',...
        'algorithm','trust-region-reflective','LargeScale','on',...
        'MaxFunEvals',10^13,'MaxIter',5000,...
        'display','notify',...
        'FunValCheck','on','DerivativeCheck','off');
    % DerivativeCheck can be set to on when you want to make sure that your
    % analytically calculated derivative agree with the derivatives
    % calculated using finite differences
    
    % Run this part in a try loop because occasionally it crashes with some Matlab versions
    if ignoreerrors
        try
            result = fmincon(@(parameters) calculateerrorMJ3D(parameters,time,v,tv,time(2)-time(1)),initialparameters,[],[],[],[],lb,ub,[],options);
            [epsilon,~,~,fitresult] = calculateerrorMJ3D(result,time,v,tv,time(2)-time(1));
            
            
            if ~isreal(result(1))
                error('Found an imaginary value');
            end
            
            if epsilon < bestError
                bestError = epsilon;
                bestParameters = result;
                bestVelocity = fitresult;
            end
            
        catch exception
            % occasionally there are errors for reasons that are not clear, let
            % us just ignore them [seem to have gone away now]
            fprintf(['Got a strange error: ' exception.message ' in file ' exception.stack(1).name ' on line ' num2str(exception.stack(1).line) ', ignoring\n']);
            count = count-1;
        end
        count = count+1;
    else
        result = fmincon(@(parameters) calculateerrorMJ3D(parameters,time,v,tv,time(2)-time(1)),initialparameters,[],[],[],[],lb,ub,[],options);
        [epsilon,~,~,fitresult] = calculateerrorMJ3D(result,time,v,tv,time(2)-time(1));
        
        
        if ~isreal(result(1))
            error('Found an imaginary value');
        end
        
        if epsilon < bestError
            bestError = epsilon;
            bestParameters = result;
            bestVelocity = fitresult;
        end
        count = count+1;
    end
end

% Sort the parameters according to t0
t0 = bestParameters(1:pps:end-4);
D  = bestParameters(2:pps:end-3);
Ax = bestParameters(3:pps:end-2);
Ay = bestParameters(4:pps:end-1);
Az = bestParameters(5:pps:end);

[~,order] = sort(t0);
t0 = t0(order);
D = D(order);
Ax = Ax(order);
Ay = Ay(order);
Az = Az(order);

bestParameters(1:pps:end-4) = t0;
bestParameters(2:pps:end-3) = D;
bestParameters(3:pps:end-2) = Ax;
bestParameters(4:pps:end-1) = Ay;
bestParameters(5:pps:end)   = Az;