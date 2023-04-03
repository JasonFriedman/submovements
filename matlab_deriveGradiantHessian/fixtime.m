% FIXTIME  - we only want to evaluate over the range r 
% (because the min jerk function is undefined outside this range)

function fixed = fixtime(s)

ts = strfind(s,'t');
for k=numel(ts):-1:1
    % don't add (r) if it is t0 or part of sqrt
    if (numel(s)==ts(k) || s(ts(k)+1)~='0') && ~strcmp(s(ts(k)-3:ts(k)),'sqrt')
        s = [s(1:ts(k)) '(r)' s(ts(k)+1:end)];
    end
end

fixed = s;