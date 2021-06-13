% PLOTPOSITION - plot the position data
% 
% plotposition(position,time,plottype)
%
% plottype = 1 (default) → x vs y
% plottype = 2           → time vs x/y

function plotposition(position,time,plottype)

if nargin<3
    plottype=1;
end

figure;
cols = ceil(sqrt(numel(position)));
rows = ceil(numel(position) / cols);

for k=1:numel(position)
        subplot(rows,cols,k);
        if plottype==1
            plot(position{k}(:,1),position{k}(:,2));
            axis equal
        elseif plottype==2
            plot(time{k},position{k});
            if k==numel(position)
                legend('x','y');
            end
        else
            error('Unknown plot type');
        end
end
