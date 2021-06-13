% CALCULATERELATIVEONSET - Calculate relative onset time of the submovements
%
% relativeonset = calculateRelativeOnset(parameters) 
%
% Values > 100% indicate that the 2nd submovement started after the 
% first one finished
%
% parameters should be an N * p matrix
% with N repetitions of the task, and p = 4 * submovements

function relativeOnsets = calculateRelativeOnset(parameters) 

if size(parameters,1)>1
    for k=1:size(parameters,1)
        relativeOnsets(k,:) = calculateRelativeOnset(parameters(k,:));
    end
    return
end

numsubmovements = numel(parameters)/4;

relativeOnsets = [];

for k=2:numsubmovements
   t01 = parameters((k-2)*4+1);
   D1  = parameters((k-2)*4+2);
   t02 = parameters((k-1)*4+1);
   D2  = parameters((k-1)*4+2);
   
   relativeOnsets(k-1) = (t02 - t01) / D1 * 100;
end