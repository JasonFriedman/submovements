% PLOTPOSITION - plot the position data (in 2D)
% 
% plotposition(position,time,plottype)
%
% position should be a cell array of N x 2 matrices with the x and y positions for each time vector
% time should be a cell array of N x 1 vectors with the corresponding time (in seconds) for each position vector
% plottype = 1 (default) → x vs y
% plottype = 2           → time vs x/y
%

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
