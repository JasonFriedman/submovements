function constraints = resolveFittingConstraints(overrides)
% RESOLVEFITTINGCONSTRAINTS - merge default fitting constraints with user overrides

defaults.minOnsetSpacing = 0.167;
defaults.minDuration = 0.167;
defaults.maxDuration = 1.0;
defaults.minUpperBoundTime = 0.1;
defaults.numRestarts = 20;
defaults.maxFunEvals = 10^13;
defaults.maxIter = 5000;

if nargin<1 || isempty(overrides)
    constraints = defaults;
    return
end

if ~isstruct(overrides)
    error('fittingConstraints must be a struct');
end

constraints = defaults;
fields = fieldnames(overrides);
for k=1:numel(fields)
    fieldName = fields{k};
    if ~isfield(defaults,fieldName)
        error('Unknown fitting constraint: %s',fieldName);
    end
    constraints.(fieldName) = overrides.(fieldName);
end

if constraints.minOnsetSpacing <= 0
    error('minOnsetSpacing must be > 0');
end
if constraints.minDuration <= 0
    error('minDuration must be > 0');
end
if constraints.maxDuration <= 0
    error('maxDuration must be > 0');
end
if constraints.maxDuration < constraints.minDuration
    error('maxDuration must be >= minDuration');
end
if constraints.minUpperBoundTime <= 0
    error('minUpperBoundTime must be > 0');
end
if constraints.numRestarts < 1 || mod(constraints.numRestarts,1)~=0
    error('numRestarts must be a positive integer');
end
if constraints.maxFunEvals < 1
    error('maxFunEvals must be >= 1');
end
if constraints.maxIter < 1 || mod(constraints.maxIter,1)~=0
    error('maxIter must be a positive integer');
end
