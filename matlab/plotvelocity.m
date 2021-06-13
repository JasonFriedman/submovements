% PLOTVELOCITY - plot the velocity data
% 
% plotvelocity(velocity,time,plottype)
%
% plottype = 1 (default) → time vs v_x/v_y
% plottype = 2           → time vs tangential velocity

function plotvelocity(velocity,time,plottype)

if nargin<3
    plottype=1;
end

figure;
cols = ceil(sqrt(numel(velocity)));
rows = ceil(numel(velocity) / cols);

for k=1:numel(velocity)
        subplot(rows,cols,k);
        if plottype==1
            plot(time{k},velocity{k});
            if k==numel(velocity)
                legend('v_x','v_y');
            end
        elseif plottype==2
            tangvel = sqrt(sum(velocity{k}.^2,2));
            plot(time{k},tangvel);
        else
            error('Unknown plot type');
        end
end
