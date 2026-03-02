% PLOTSUBMOVEMENTS2D - plot 2D submovements after decomposition
%
% plotSubmovements2D(parameters,t,plottype,x0,y0)
%
% The parameters should be in sets of 4 for each submovement:
% [t0 D Ax Ay]
%
%
% plottype:
% 1 = time vs submovement velocity + sum velocity (default)
% 2 = time vs submovement velocity
% 3 = time vs submovement position - extra parameters (x0,y0) specifies the start
%     position of the first submovement (the other submovements are assumed
%     to start where the previous submovement ended)
% 4 = same as 3, but without the sum
% 5 = submovement position x vs y + sum

function plotSubmovements2D(parameters,t,plottype,x0,y0)

if nargin<3 || isempty(plottype)
    plottype=1;
end

if nargin<4 || isempty(x0)
    x0 = 0;
end

if nargin<5 || isempty(y0)
    y0 = 0;
end

if nargin<2 || isempty(t)
    t = [];
end

data = prepareSubmovementKinematics(parameters,2,t,[x0 y0]);
t = data.t;
t0 = data.t0;
D = data.D;
numSubmovements = data.numSubmovements;
vx = reshape(data.vel(:,:,1),numSubmovements,[]);
vy = reshape(data.vel(:,:,2),numSubmovements,[]);
x = reshape(data.pos(:,:,1),numSubmovements,[]);
y = reshape(data.pos(:,:,2),numSubmovements,[]);
x_sum = data.sumPos(:,1)';
y_sum = data.sumPos(:,2)';

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