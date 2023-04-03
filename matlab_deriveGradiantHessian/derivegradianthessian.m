%% This script derives the gradiant and Hessian
% The printed values are then copy / pasted into minimumJerkVelocity2D / 3D

syms t t0 D Ax Ay Az

x = Ax * (6*((t-t0)/D)^5 - 15*((t-t0)/D)^4 + 10*((t-t0)/D)^3);
y = Ay * (6*((t-t0)/D)^5 - 15*((t-t0)/D)^4 + 10*((t-t0)/D)^3);
z = Az * (6*((t-t0)/D)^5 - 15*((t-t0)/D)^4 + 10*((t-t0)/D)^3);

xdot = diff(x,t);
ydot = diff(y,t);
zdot = diff(z,t);

%% 2D
tangvel = simplify(sqrt(xdot^2 + ydot^2));
varnames = {'t0','D','Ax','Ay'};
toderive = {'xdot','ydot','tangvel'};
outputnames = {'Jx','Jy','J'};

%% gradiant (partial derivatives) - 2D
for outputVar = 1:numel(toderive)
    for varname = 1:numel(varnames)
        eval(['s = diff(' toderive{outputVar} ',' varnames{varname} ');'])
        fprintf('%s(%d,r) = %s;\n',outputnames{outputVar},varname,fixtime(cleanExpression(s)));
    end
end

%% gradiant (partial derivitives) - 2D

outputnames = {'Hx','Hy','H'};

for outputVar = 1:numel(outputnames)
    for varname1 = 1:numel(varnames)
        for varname2 = 1:numel(varnames)
            eval(['s = simplify(diff(diff(' toderive{outputVar} ',' varnames{varname1} '),' varnames{varname2} '));'])
            fprintf('%s(%d,%d,r) = %s;\n',outputnames{outputVar},varname1,varname2,fixtime(cleanExpression(s)));
        end
    end
end

%% 3D
tangvel = simplify(sqrt(xdot^2 + ydot^2 + zdot^2));
varnames = {'t0','D','Ax','Ay','Az'};
toderive = {'xdot','ydot','zdot','tangvel'};
outputnames = {'Jx','Jy','Jz','J'};
%% gradiant (partial derivitives) - 3D
for outputVar = 1:numel(toderive)
    for varname = 1:numel(varnames)
        eval(['s = diff(' toderive{outputVar} ',' varnames{varname} ');'])
        fprintf('%s(%d,r) = %s;\n',outputnames{outputVar},varname,fixtime(cleanExpression(s)));
    end
end

%% Hessian - 3D

outputnames = {'Hx','Hy','Hz','H'};

for outputVar = 1:numel(outputnames)
    for varname1 = 1:numel(varnames)
        for varname2 = 1:numel(varnames)
            eval(['s = simplify(diff(diff(' toderive{outputVar} ',' varnames{varname1} '),' varnames{varname2} '));'])
            fprintf('%s(%d,%d,r) = %s;\n',outputnames{outputVar},varname1,varname2,fixtime(cleanExpression(s)));
        end
    end
end



