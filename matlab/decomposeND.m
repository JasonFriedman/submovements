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

bestError = inf;
ignoreerrors = true;

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

count=1;
while count<=constraints.numRestarts
    for i=1:numsubmovements
        initialparameters(1,i*parametersPerSubmovement-(parametersPerSubmovement-1):i*parametersPerSubmovement) = ...
            lb_0 + (ub_0-lb_0) .* rand(1,parametersPerSubmovement);
    end

    options = optimset('GradObj','on','Hessian','on',...
        'algorithm','trust-region-reflective','LargeScale','on',...
        'MaxFunEvals',constraints.maxFunEvals,'MaxIter',constraints.maxIter,...
        'display','notify',...
        'FunValCheck','on','DerivativeCheck','off');

    if ignoreerrors
        try
            result = fmincon(errorFunction,initialparameters,[],[],[],[],lb,ub,[],options);
            [epsilon,~,~,fitresult] = errorFunction(result);

            if ~isreal(result(1))
                error('Found an imaginary value');
            end

            if epsilon < bestError
                bestError = epsilon;
                bestParameters = result;
                bestVelocity = fitresult;
            end

        catch exception
            fprintf(['Got a strange error: ' exception.message ' in file ' exception.stack(1).name ' on line ' num2str(exception.stack(1).line) ', ignoring\n']);
            count = count-1;
        end
        count = count+1;
    else
        result = fmincon(errorFunction,initialparameters,[],[],[],[],lb,ub,[],options);
        [epsilon,~,~,fitresult] = errorFunction(result);

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

[~,order] = sort(bestParameters(1:parametersPerSubmovement:end-parametersPerSubmovement+1));
for parameterIndex=1:parametersPerSubmovement
    values = bestParameters(parameterIndex:parametersPerSubmovement:end-parametersPerSubmovement+parameterIndex);
    bestParameters(parameterIndex:parametersPerSubmovement:end-parametersPerSubmovement+parameterIndex) = values(order);
end
