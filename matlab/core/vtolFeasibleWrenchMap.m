function map = vtolFeasibleWrenchMap(p, architecture, gridSpec)
%VTOLFEASIBLEWRENCHMAP Constrained nonlinear authority maps about hover.
%
%   map = vtolFeasibleWrenchMap(p,architecture,gridSpec) allocates each
%   requested wrench while holding the remaining axes at hover values.  It
%   returns requested-versus-achieved Fx/My and Fy/My grids, the resulting
%   feasible envelopes, maximum angular acceleration per axis, and scaled
%   Jacobian metrics over representative gimbal grids.
%
%   A grid point is feasible only if the weighted allocation residual is
%   below gridSpec.scaledResidualTolerance and Fz remains within
%   gridSpec.fzTolerance of -m*g.  This is a numerical constrained map, not
%   a hover-only linear approximation.

if nargin < 3 || isempty(gridSpec)
    gridSpec = struct();
end
key = vtolArchitectureKey(architecture);
spec = localGridSpec(p, gridSpec);
[uTrim, Wtrim, trimInfo] = vtolHoverTrim(p, key);
Whold = [0; 0; -p.m*p.g; 0; 0; 0];

allocationOptions = struct();
allocationOptions.solver = spec.solver;
allocationOptions.maxIterations = spec.maxIterations;
allocationOptions.enforceRateBounds = false;
allocationOptions.feasibilityTolerance = spec.scaledResidualTolerance;

fxRequested = linspace(-spec.FxLimit, spec.FxLimit, spec.nFx);
fyRequested = linspace(-spec.FyLimit, spec.FyLimit, spec.nFy);
myRequested = linspace(-spec.MyLimit, spec.MyLimit, spec.nMy);

map = struct();
map.architecture = key;
map.parameters = p;
map.trim = struct('u',uTrim,'W',Wtrim,'info',trimInfo);
map.specification = spec;
map.fxMy = localPairGrid(fxRequested, myRequested, 1, 5, ...
    Whold, uTrim, p, key, allocationOptions, spec);
map.fyMy = localPairGrid(fyRequested, myRequested, 2, 5, ...
    Whold, uTrim, p, key, allocationOptions, spec);
map.axisAuthority = localAxisAuthority(Whold, uTrim, p, key, ...
    allocationOptions, spec);
map.jacobian = localJacobianGrid(uTrim, p, key, spec);
map.scalingNote = ['Jscaled = diag(1./p.scaling.wrench)*J*diag(half actuator travel). ' ...
    'Condition numbers in this map are only for Jscaled.'];
end

function grid = localPairGrid(primaryRequested, secondaryRequested, ...
    primaryIndex, secondaryIndex, Whold, uTrim, p, key, allocationOptions, spec)

nPrimary = numel(primaryRequested);
nSecondary = numel(secondaryRequested);
nActuator = numel(uTrim);

grid = struct();
grid.primaryIndex = primaryIndex;
grid.secondaryIndex = secondaryIndex;
grid.primaryRequested = primaryRequested;
grid.secondaryRequested = secondaryRequested;
grid.actual = zeros(6,nSecondary,nPrimary);
grid.command = zeros(nActuator,nSecondary,nPrimary);
grid.residual = zeros(6,nSecondary,nPrimary);
grid.scaledResidualNorm = zeros(nSecondary,nPrimary);
grid.uMax = zeros(nSecondary,nPrimary);
grid.utilization = zeros(nActuator,nSecondary,nPrimary);
grid.activeConstraintCount = zeros(nSecondary,nPrimary);
grid.feasible = false(nSecondary,nPrimary);

for i = 1:nPrimary
    for j = 1:nSecondary
        Wdes = Whold;
        Wdes(primaryIndex) = primaryRequested(i);
        Wdes(secondaryIndex) = secondaryRequested(j);
        result = vtolConstrainedAllocate(Wdes, uTrim, p, key, allocationOptions);

        grid.actual(:,j,i) = result.W_actual;
        grid.command(:,j,i) = result.u_cmd;
        grid.residual(:,j,i) = result.residual;
        grid.scaledResidualNorm(j,i) = result.scaledResidualNorm;
        grid.uMax(j,i) = result.u_max;
        grid.utilization(:,j,i) = result.utilization;
        grid.activeConstraintCount(j,i) = sum(result.activeConstraints.any);
        grid.feasible(j,i) = result.scaledResidualNorm <= ...
            spec.scaledResidualTolerance && ...
            abs(result.W_actual(3) + p.m*p.g) <= spec.fzTolerance;
    end
end

grid.envelope = localEnvelopes(grid);
end

function envelope = localEnvelopes(grid)
nPrimary = numel(grid.primaryRequested);
nSecondary = numel(grid.secondaryRequested);
secondaryMinimum = nan(1,nPrimary);
secondaryMaximum = nan(1,nPrimary);
primaryMinimum = nan(1,nSecondary);
primaryMaximum = nan(1,nSecondary);

for i = 1:nPrimary
    feasibleRows = find(grid.feasible(:,i));
    if ~isempty(feasibleRows)
        value = reshape(grid.actual(grid.secondaryIndex,feasibleRows,i),[],1);
        secondaryMinimum(i) = min(value);
        secondaryMaximum(i) = max(value);
    end
end
for j = 1:nSecondary
    feasibleColumns = find(grid.feasible(j,:));
    if ~isempty(feasibleColumns)
        value = reshape(grid.actual(grid.primaryIndex,j,feasibleColumns),[],1);
        primaryMinimum(j) = min(value);
        primaryMaximum(j) = max(value);
    end
end

envelope = struct();
envelope.secondaryMinimumVsPrimaryRequest = secondaryMinimum;
envelope.secondaryMaximumVsPrimaryRequest = secondaryMaximum;
envelope.primaryMinimumVsSecondaryRequest = primaryMinimum;
envelope.primaryMaximumVsSecondaryRequest = primaryMaximum;
end

function authority = localAxisAuthority(Whold, uTrim, p, key, allocationOptions, spec)
axisNames = {'Mx','My','Mz'};
authority = struct();
authority.names = axisNames;
authority.requested = cell(1,3);
authority.actual = cell(1,3);
authority.feasible = cell(1,3);
authority.minimum = nan(3,1);
authority.maximum = nan(3,1);
authority.angularAccelerationMinimum = nan(3,1);
authority.angularAccelerationMaximum = nan(3,1);
authority.angularAccelerationMagnitude = nan(3,1);

for k = 1:3
    requested = linspace(-spec.MomentLimits(k), spec.MomentLimits(k), spec.nAxis);
    actual = nan(size(requested));
    feasible = false(size(requested));
    for j = 1:numel(requested)
        Wdes = Whold;
        Wdes(3+k) = requested(j);
        result = vtolConstrainedAllocate(Wdes, uTrim, p, key, allocationOptions);
        actual(j) = result.W_actual(3+k);
        feasible(j) = result.scaledResidualNorm <= spec.scaledResidualTolerance && ...
            abs(result.W_actual(3)+p.m*p.g) <= spec.fzTolerance;
    end
    authority.requested{k} = requested;
    authority.actual{k} = actual;
    authority.feasible{k} = feasible;
    if any(feasible)
        authority.minimum(k) = min(actual(feasible));
        authority.maximum(k) = max(actual(feasible));
        authority.angularAccelerationMinimum(k) = authority.minimum(k)/p.I(k,k);
        authority.angularAccelerationMaximum(k) = authority.maximum(k)/p.I(k,k);
        authority.angularAccelerationMagnitude(k) = max(abs([ ...
            authority.angularAccelerationMinimum(k), ...
            authority.angularAccelerationMaximum(k)]));
    end
end
end

function jacobian = localJacobianGrid(uTrim, p, key, spec)
alphaGrid = linspace(-p.main.alphaMax,p.main.alphaMax,spec.nAngle);
betaGrid = linspace(-p.main.betaMax,p.main.betaMax,spec.nAngle);
modeNames = {'commonGimbal','differentialBeta'};

jacobian = struct();
jacobian.alphaGrid = alphaGrid;
jacobian.betaGrid = betaGrid;
jacobian.mode = struct( ...
    'name',{}, ...
    'rank',{}, ...
    'conditionNumber',{}, ...
    'singularValues',{});

for modeIndex = 1:numel(modeNames)
    mode = struct();
    mode.name = modeNames{modeIndex};
    mode.rank = zeros(numel(betaGrid),numel(alphaGrid));
    mode.conditionNumber = nan(numel(betaGrid),numel(alphaGrid));
    mode.singularValues = zeros(6,numel(betaGrid),numel(alphaGrid));
    for ia = 1:numel(alphaGrid)
        for ib = 1:numel(betaGrid)
            u = uTrim;
            u([2 5]) = alphaGrid(ia);
            if modeIndex == 1
                u([3 6]) = betaGrid(ib);
            else
                u(3) = betaGrid(ib);
                u(6) = -betaGrid(ib);
            end
            J = vtolWrenchJacobian(u,p,key,'analytic');
            metrics = vtolScaleJacobian(J,p,key);
            mode.rank(ib,ia) = metrics.rank;
            mode.conditionNumber(ib,ia) = metrics.conditionNumber;
            mode.singularValues(:,ib,ia) = metrics.singularValues;
        end
    end
    jacobian.mode(modeIndex) = mode;
end
jacobian.note = ['commonGimbal applies alpha1=alpha2 and beta1=beta2; ' ...
    'differentialBeta applies alpha1=alpha2 and beta1=-beta2.'];
end

function spec = localGridSpec(p, supplied)
spec = struct();
spec.nFx = localField(supplied,'nFx',p.map.nFx);
spec.nFy = localField(supplied,'nFy',p.map.nFy);
spec.nMy = localField(supplied,'nMy',p.map.nMy);
spec.nAxis = localField(supplied,'nAxis',p.map.nAxis);
spec.nAngle = localField(supplied,'nAngle',p.map.nAngle);
spec.FxLimit = localField(supplied,'FxLimit',p.map.FxLimit);
spec.FyLimit = localField(supplied,'FyLimit',p.map.FyLimit);
spec.MyLimit = localField(supplied,'MyLimit',p.map.MyLimit);
spec.MomentLimits = localField(supplied,'MomentLimits',p.map.MomentLimits);
spec.fzTolerance = localField(supplied,'fzTolerance',p.map.fzTolerance);
spec.scaledResidualTolerance = localField(supplied, ...
    'scaledResidualTolerance',p.map.scaledResidualTolerance);
spec.solver = localField(supplied,'solver',p.map.solver);
spec.maxIterations = localField(supplied,'maxIterations',p.map.maxIterations);

if any([spec.nFx spec.nFy spec.nMy spec.nAxis spec.nAngle] < 2)
    error('vtolFeasibleWrenchMap:GridSize', ...
        'Every requested grid dimension must contain at least two samples.');
end
end

function value = localField(s, name, defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
