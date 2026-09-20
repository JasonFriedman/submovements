function [bestError,bestParameters,bestVelocity] = decomposeND(time,numsubmovements,lb_0,ub_0,parametersPerSubmovement,errorFunction,fittingConstraints)
% DECOMPOSEND - shared decomposition engine for 2D/3D minimum jerk models

if nargin<7
    fittingConstraints = struct();
end

constraints = resolveFittingConstraints(fittingConstraints);

if numel(time)==0
    bestError = NaN;
    bestParameters = NaN(1,numsubmovements*parametersPerSubmovement);
    bestVelocity = NaN;
    return
end

ignoreerrors = true;

if numsubmovements > 1
    nonlcon = @(parameters) onsetSpacingConstraint(parameters,constraints.minOnsetSpacing,parametersPerSubmovement,numsubmovements);
else
    nonlcon = [];
end

if any(lb_0>ub_0)
    error('Lower bounds exceed upper bound - infeasible');
end

toignore = false;
for i=1:numsubmovements
    thislb_0 = lb_0;
    thislb_0(1) = (i-1)*constraints.minOnsetSpacing;
    if thislb_0(1) > ub_0(1)
        fprintf(['The submovements are assumed to be spaced by at least 167 ms, this movement is not long enough for ' num2str(i) ' submovements so will be set to NaN\n']);
        toignore = true;
        break;
    end
    lb(i*parametersPerSubmovement-(parametersPerSubmovement-1):i*parametersPerSubmovement) = thislb_0;
    ub(i*parametersPerSubmovement-(parametersPerSubmovement-1):i*parametersPerSubmovement) = ub_0;
end

if toignore
    bestError = NaN;
    bestParameters = NaN;
    bestVelocity = NaN;
    return;
end

% options only depend on whether nonlcon is present, not on the restart
% index, so build them once rather than on every restart iteration.
if isempty(nonlcon)
    options = optimset('GradObj','on','Hessian','on',...
        'algorithm','trust-region-reflective','LargeScale','on',...
        'MaxFunEvals',constraints.maxFunEvals,'MaxIter',constraints.maxIter,...
        'display','notify',...
        'FunValCheck','on','DerivativeCheck','off');
else
    % trust-region-reflective does not support nonlinear constraints.
    options = optimset('GradObj','on','Hessian','off',...
        'algorithm','interior-point','LargeScale','on',...
        'MaxFunEvals',constraints.maxFunEvals,'MaxIter',constraints.maxIter,...
        'display','notify',...
        'FunValCheck','on','DerivativeCheck','off');
end

% Restarts are independent, so they run in a parfor. If Parallel Computing
% Toolbox is not installed/licensed, parfor automatically falls back to
% running the loop body serially (like a normal for), so no separate
% code path is needed.
numRestarts = constraints.numRestarts;
restartEpsilon = nan(1,numRestarts);
restartParameters = cell(1,numRestarts);
restartVelocity = cell(1,numRestarts);

parfor count=1:numRestarts
    initialparameters = zeros(1,numsubmovements*parametersPerSubmovement);
    for i=1:numsubmovements
        cols = i*parametersPerSubmovement-(parametersPerSubmovement-1):i*parametersPerSubmovement;
        initialparameters(cols) = lb_0 + (ub_0-lb_0) .* rand(1,parametersPerSubmovement);
    end

    if ignoreerrors
        try
            result = fmincon(errorFunction,initialparameters,[],[],[],[],lb,ub,nonlcon,options);
            [epsilon,~,~,fitresult] = errorFunction(result);

            if ~isreal(result(1))
                error('Found an imaginary value');
            end

            restartEpsilon(count) = epsilon;
            restartParameters{count} = result;
            restartVelocity{count} = fitresult;
        catch exception
            fprintf(['Got a strange error: ' exception.message ' in file ' exception.stack(1).name ' on line ' num2str(exception.stack(1).line) ', ignoring\n']);
        end
    else
        result = fmincon(errorFunction,initialparameters,[],[],[],[],lb,ub,nonlcon,options);
        [epsilon,~,~,fitresult] = errorFunction(result);

        if ~isreal(result(1))
            error('Found an imaginary value');
        end

        restartEpsilon(count) = epsilon;
        restartParameters{count} = result;
        restartVelocity{count} = fitresult;
    end
end

[bestError,bestIdx] = min(restartEpsilon);
if isnan(bestError)
    error('All %d restarts failed to produce a valid fit (numsubmovements=%d)',numRestarts,numsubmovements);
end
bestParameters = restartParameters{bestIdx};
bestVelocity = restartVelocity{bestIdx};

[~,order] = sort(bestParameters(1:parametersPerSubmovement:end-parametersPerSubmovement+1));
for parameterIndex=1:parametersPerSubmovement
    values = bestParameters(parameterIndex:parametersPerSubmovement:end-parametersPerSubmovement+parameterIndex);
    bestParameters(parameterIndex:parametersPerSubmovement:end-parametersPerSubmovement+parameterIndex) = values(order);
end

function [c,ceq] = onsetSpacingConstraint(parameters,minOnsetSpacing,parametersPerSubmovement,numsubmovements)
% Enforce pairwise onset spacing: t0(i+1)-t0(i) >= minOnsetSpacing.
t0s = parameters(1:parametersPerSubmovement:numsubmovements*parametersPerSubmovement);
if numel(t0s) <= 1
    c = [];
else
    c = minOnsetSpacing - diff(t0s);
end
ceq = [];
