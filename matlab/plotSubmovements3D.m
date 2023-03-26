% PLOTSUBMOVEMENTS3D - plot 3D submovements after decomposition
%
% plotSubmovements3D(parameters,t,plottype,x0,y0)
%
% The parameters should be in sets of 5 for each submovement:
% [t0 D Ax Ay Az]
%
%
% plottype:
% 1 = time vs submovement velocity + sum velocity (default)
% 2 = time vs submovement velocity
% 3 = time vs submovement position - extra parameters (x0,y0,z0) specifies the start
%     position of the first submovement (the other submovements are assumed
%     to start where the previous submovement ended)
% 4 = same as 3, but without the sum
% 5 = submovement position x vs y + sum

function plotSubmovements3D(parameters,t,plottype,x0,y0)

if mod(numel(parameters),4)~=0
    error('The parameters vector must have a length that is a multiple of 4');
end

if nargin<3 || isempty(plottype)
    plottype=1;
end

if nargin<4 || isempty(x0)
    x0 = 0;
end

if nargin<5 || isempty(y0)
    y0 = 0;
end

numSubmovements = numel(parameters)/4;
t0 = parameters(1:5:end-4);
D =  parameters(2:5:end-3);
Ax = parameters(3:5:end-2);
Ay = parameters(4:5:end-1);
Az = parameters(5:5:end);

% sort according to t0
[~,order] = sort(t0);
t0 = t0(order);
D = D(order);
Ax = Ax(order);
Ay = Ay(order);
Az = Az(order);

x0(2:numSubmovements) = x0 + cumsum(Ax(1:end-1));
y0(2:numSubmovements) = y0 + cumsum(Ay(1:end-1));
z0(2:numSubmovements) = z0 + cumsum(Az(1:end-1));

tf = t0+D;

if nargin<2 || isempty(t)
    t = linspace(min(t0),max(tf),100);
end

for s = 1:numSubmovements    
    [vx(s,:),vy(s,:),vz(s,:)] = minimumJerkVelocity2D(t0(s),D(s),Ax(s),Ay(s),t);
    [x(s,:),y(s,:),z(s,:)] = minimumJerkPosition2D(t0(s),D(s),Ax(s),Ay(s),x0(s),y0(s),t);
end

x_relative(1,:) = x(1,:); y_relative(1,:) = y(1,:); z_relative(1,:) = z(1,:);
for s=2:numSubmovements
    x_relative(s,:) = x(s,:) - x0(s);
    y_relative(s,:) = y(s,:) - y0(s);
    z_relative(s,:) = z(s,:) - z0(s);
end
x_sum = sum(x_relative);
y_sum = sum(y_relative);
z_sum = sum(z_relative);

if any(plottype==1:2)
    h(1:2:numSubmovements*2-1) = plot(t,vx,'b');
    hold on;
    h(2:2:numSubmovements*2) = plot(t,vy,'r');
    legend(h(1:2),'Submovements v_x','Submovements v_y');
    xlabel('time');
    ylabel('velocity');
end

if any(plottype==3:4)
    for s=1:size(x,1)
        r = find(t>=t0(s) & t<=t0(s)+D(s));
        h(1) = plot(t(r),x(s,r),'b');
        hold on;
        h(2) = plot(t(r),y(s,r),'r');
    end
    legend(h(1:2),'Submovements x','Submovements y');
    xlabel('time');
    ylabel('position');
end

if plottype==1
    hh(1) = plot(t,sum(vx),'k--','LineWidth',2);
    hh(2) = plot(t,sum(vy),'m--','LineWidth',2);
    legend([h(1:2) hh(1:2)],'Submovements v_x','Submovements v_y','Sum movements v_x','Sum movements v_y');
end

if plottype==3
    hh(1) = plot(t,x_sum,'k--','LineWidth',2);
    hh(2) = plot(t,y_sum,'m--','LineWidth',2);
    legend([h(1:2) hh(1:2)],'Submovements x','Submovements y','Sum submovements x','Sum submovements y');
end

if plottype==5
    for s=1:size(x,1)
        r = find(t>=t0(s) & t<=t0(s)+D(s));
        h(1) = plot(x(s,r),y(s,r),'b');
        hold on;
    end
    hh(1) = plot(x_sum,y_sum,'k--','LineWidth',2);
    xlabel('x');
    ylabel('y');
    legend([h(1) hh(1)],'Submovements','Sum submovements');
end