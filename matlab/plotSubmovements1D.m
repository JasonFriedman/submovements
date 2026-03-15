% PLOTSUBMOVEMENTS1D - plot 1D submovements after decomposition
%
% plotSubmovements1D(parameters,t,plottype,x0)
%
% The parameters should be in sets of 3 for each submovement:
% [t0 D A]
%
%
% plottype:
% 1 = time vs submovement velocity + sum velocity (default)
% 2 = time vs submovement velocity
% 3 = time vs submovement position - extra parameter x0 specifies the start
%     position of the first submovement (the other submovements are assumed
%     to start where the previous submovement ended)
% 4 = same as 3, but without the sum

function plotSubmovements1D(parameters,t,plottype,x0)

if nargin<3 || isempty(plottype)
    plottype=1;
end

if nargin<4 || isempty(x0)
    x0 = 0;
end

if nargin<2 || isempty(t)
    t = [];
end

if mod(numel(parameters),3)~=0
    error('The parameters vector must have a length that is a multiple of 3');
end

numSubmovements = numel(parameters)/3;
parameterMatrix = reshape(parameters,3,numSubmovements)';
t0 = parameterMatrix(:,1)';
D = parameterMatrix(:,2)';
A = parameterMatrix(:,3)';

[~,order] = sort(t0);
t0 = t0(order);
D = D(order);
A = A(order);

starts = zeros(1,numSubmovements);
starts(1) = x0;
if numSubmovements>1
    starts(2:numSubmovements) = x0 + cumsum(A(1:end-1));
end

tf = t0 + D;
if isempty(t)
    t = linspace(min(t0),max(tf),100);
end

v = zeros(numSubmovements,numel(t));
x = zeros(numSubmovements,numel(t));
for s=1:numSubmovements
    v(s,:) = minimumJerkVelocity1D(t0(s),D(s),A(s),t);
    x(s,:) = minimumJerkPosition1D(t0(s),D(s),A(s),starts(s),t);
end

xRelative = x;
for s=2:numSubmovements
    xRelative(s,:) = x(s,:) - starts(s);
end
x_sum = sum(xRelative,1)';

if any(plottype==1:2)
    h = plot(t,v,'b');
    xlabel('time');
    ylabel('velocity');
    legend(h(1),'Submovements velocity');
end

if any(plottype==3:4)
    holdState = ishold;
    for s=1:size(x,1)
        r = find(t>=t0(s) & t<=t0(s)+D(s));
        h = plot(t(r),x(s,r),'b');
        hold on;
    end
    if ~holdState && plottype~=3
        hold off;
    end
    xlabel('time');
    ylabel('position');
    legend(h(1),'Submovements x');
end

if plottype==1
    hold on;
    hh = plot(t,sum(v),'k--','LineWidth',2);
    legend([h(1) hh(1)],'Submovements velocity','Sum movements velocity');
end

if plottype==3
    hold on;
    hh = plot(t,x_sum,'k--','LineWidth',2);
    legend([h(1) hh(1)],'Submovements x','Sum submovements x');
end