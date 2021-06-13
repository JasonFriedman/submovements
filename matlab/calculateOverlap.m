% CALCULATEOVERLAP - Calculate temporal overlap between submovements
%
% overlap = calculateOverlap(parameters) 
%
% parameters should be an N * p matrix
% with N repetitions of the task, and p = 4 * submovements

function overlaps = calculateOverlap(parameters)

if size(parameters,1)>1
    for k=1:size(parameters,1)
        overlaps(k,:) = calculateOverlap(parameters(k,:));
    end
    return
end

numsubmovements = numel(parameters)/4;

overlaps = [];

for k=2:numsubmovements
   t01 = parameters((k-2)*4+1);
   D1  = parameters((k-2)*4+2);
   t02 = parameters((k-1)*4+1);
   D2  = parameters((k-1)*4+2);

   if t02 > t01 + D1
       overlaps(k-1) = 0;
   else
       overlaps(k-1) = min([(t01 + D1 - t02) / (t02 + D2 - t01) * 100, 100]);
   end
end