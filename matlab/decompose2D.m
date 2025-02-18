function [bestError,bestParameters,bestVelocity] = decompose2D(time,vel,numsubmovements,xrng,yrng,criteria)
% DECOMPOSE2D - decompose two dimensional movement into submovements using the velocity profiles
%
% [best,bestParameters,bestVelocity] = decompose(time,vel,numsubmovements,xrng,yrng,criteria)
%
% vel should be a N x 2 matrix, with the x and y velocities
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
% min(t0) = 0.167 * submovement number
%
% criteria - stop if the bestError is less than the criteria (only relevant
% when numsubmovements is a vector of multiple values). This can save time
% by not checking for a higher number of submovements
%
%
% bestError the best (lowest) value of the error function
%
% bestParameters contains the function parameters corresponding to the best values
% [t0 D Ax Ay]. If there are multiple submovements, it will be have a
% length of 4*numsubmovements
%
% bestVelocity is the velocity profile coresponding to the best values

% Jason Friedman, 2021
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

if nargin<6 || isempty(criteria)
    criteria = -inf;
end

if size(time,2)>1
    error('time must be a N*1 vector');
end

if size(vel,2) >2 || size(vel,1)==2
    error('vel must be an N*2 matrix');
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
    
    for k=1:numel(numsubmovements)
        [bestError(k),bestParameters{k},bestVelocity{k}] = decompose2D(time,vel,numsubmovements(k),xrng,yrng);
        if bestError(k)<criteria
            return
        end
    end
    return
end

if numel(time)==0
    bestError = NaN;
    bestParameters = NaN(1,numsubmovements*4);
    bestVelocity = NaN;
    return
end
    

bestError = inf;
ignoreerrors = true;

% parameters are T0, D, Ax Ay
% ranges are
% 0 <= T0 <= finaltime - 0.167
% 0.167 <= D <= finaltime
% xrng(1) <= Ax <= xrng(2)
% yrng(1) <= Ay <= yrng(2)
lb_0 = [0                           0.167     xrng(1) yrng(1)];
ub_0 = [max([time(end)-0.167 0.1])  1.0       xrng(2) yrng(2)];
pps = 4; % parameters per submovement

if any(lb_0>ub_0)
    error('Lower bounds exceed upper bound - infeasible');
end

toignore = 0;
for i=1:numsubmovements
    thislb_0 = lb_0;
    thislb_0(1) = (i-1)*0.167;
    if thislb_0(1) > ub_0(1)
        fprintf(['The submovements are assumed to be spaced by at least 167 ms, this movement is not long enough for ' num2str(i) ' submovements so will be set to NaN\n']);
        toignore = 1;
        break;
    end
    lb(i*pps-(pps-1):i*pps) = thislb_0;
    ub(i*pps-(pps-1):i*pps) = ub_0;
end

if toignore
    bestError = NaN;
    bestParameters = NaN;
    bestVelocity = NaN;
    return;
end

% In Roher & Hogan 2006, they selected 10 random parameter selections
% Here we use 20 (increases a lot the likelihood to converge to the same
% solution on multiple runs)
count=1;
while count<=20
    for i=1:numsubmovements
        % Randomly select 20 starting positions in the legal range
        initialparameters(1,i*pps-(pps-1):i*pps) = lb_0 + (ub_0-lb_0) .* rand(1,pps);
    end
    v = vel(:,1:2);
    tv = sqrt(vel(:,1).^2 + vel(:,2).^2);
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
            result = fmincon(@(parameters) calculateerrorMJ2D(parameters,time,v,tv),initialparameters,[],[],[],[],lb,ub,[],options);
            [epsilon,~,~,fitresult] = calculateerrorMJ2D(result,time,v,tv);
            
            
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
        result = fmincon(@(parameters) calculateerrorMJ2D(parameters,time,v,tv),initialparameters,[],[],[],[],lb,ub,[],options);
        [epsilon,~,~,fitresult] = calculateerrorMJ2D(result,time,v,tv);
        
        
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
t0 = bestParameters(1:pps:end-pps+1);
D  = bestParameters(2:pps:end-pps+2);
Ax = bestParameters(3:pps:end-pps+3);
Ay = bestParameters(4:pps:end-pps+4);

[~,order] = sort(t0);
t0 = t0(order);
D = D(order);
Ax = Ax(order);
Ay = Ay(order);

bestParameters(1:pps:end-pps+1) = t0;
bestParameters(2:pps:end-pps+2) = D;
bestParameters(3:pps:end-pps+3) = Ax;
bestParameters(4:pps:end-pps+4) = Ay;