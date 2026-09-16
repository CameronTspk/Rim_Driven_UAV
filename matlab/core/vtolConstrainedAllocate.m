function result = vtolConstrainedAllocate(Wdes, uPrev, p, architecture, options)
%VTOLCONSTRAINEDALLOCATE Bounded nonlinear allocation from wrench to actuators.
%
%   result = vtolConstrainedAllocate(Wdes,uPrev,p,architecture,options)
%   minimizes
%
%     ||Sw*(W(u)-Wdes)||^2 + lambdaU*||Su*(z-zPrev)||^2,
%
%   subject to physical position bounds and, optionally, rate bounds.  z is
%   the actuator vector normalized by half actuator travel.  This makes the
%   projected Gauss-Newton fallback well conditioned enough for a first
%   Simulink allocator.  For offline studies, set options.solver='fmincon'.
%
%   Wdes is [Fx Fy Fz Mx My Mz].' in body coordinates.  The returned
%   W_actual is evaluated by the exact nonlinear vtolWrench model.

if nargin < 5 || isempty(options)
    options = struct();
end

key = vtolArchitectureKey(architecture);
[physicalLower, physicalUpper, rateBound, names] = vtolActuatorBounds(p, key);
Wdes = Wdes(:);
uPrev = min(max(uPrev(:), physicalLower), physicalUpper);

if numel(Wdes) ~= 6 || numel(uPrev) ~= numel(physicalLower)
    error('vtolConstrainedAllocate:InputSize', ...
        'Wdes must have six elements and uPrev must match the architecture.');
end

opts = localOptions(p, options, physicalLower, physicalUpper);
[lowerBound, upperBound] = localEffectiveBounds( ...
    physicalLower, physicalUpper, rateBound, uPrev, opts);

if any(lowerBound > upperBound)
    error('vtolConstrainedAllocate:InconsistentBounds', ...
        'Rate and position bounds have an empty intersection.');
end

% Normalize every decision variable by its full physically allowed half-span.
uCenter = 0.5*(physicalLower + physicalUpper);
uScale = 0.5*(physicalUpper - physicalLower);
zLower = (lowerBound-uCenter)./uScale;
zUpper = (upperBound-uCenter)./uScale;
zPrev = (uPrev-uCenter)./uScale;

if isfield(options,'uInitial') && ~isempty(options.uInitial)
    uInitial = min(max(options.uInitial(:), lowerBound), upperBound);
else
    uInitial = uPrev;
end
zInitial = (uInitial-uCenter)./uScale;

Sw = opts.Sw;
Su = opts.Su;
Qw = Sw.'*Sw;
Qu = Su.'*Su;
solver = lower(char(opts.solver));

if strcmp(solver,'auto')
    if exist('fmincon','file') == 2
        solver = 'fmincon';
    else
        solver = 'projectedgn';
    end
end

status = struct('solver',solver,'exitflag',NaN,'iterations',0, ...
    'message','','usedFallback',false);

switch solver
    case 'fmincon'
        if exist('fmincon','file') == 2
            try
                [zCommand, status] = localFmincon( ...
                    zInitial, zLower, zUpper, zPrev, uCenter, uScale, ...
                    Wdes, p, key, Qw, Qu, opts, status);
            catch ME
                [zCommand, status] = localProjectedGaussNewton( ...
                    zInitial, zLower, zUpper, zPrev, uCenter, uScale, ...
                    Wdes, p, key, Qw, Qu, opts, status);
                status.usedFallback = true;
                status.message = sprintf('fmincon failed (%s); %s', ...
                    ME.identifier, status.message);
            end
        else
            [zCommand, status] = localProjectedGaussNewton( ...
                zInitial, zLower, zUpper, zPrev, uCenter, uScale, ...
                Wdes, p, key, Qw, Qu, opts, status);
            status.usedFallback = true;
            status.message = ['fmincon unavailable; ' status.message];
        end

    case {'projectedgn','projected-gauss-newton','realtime'}
        [zCommand, status] = localProjectedGaussNewton( ...
            zInitial, zLower, zUpper, zPrev, uCenter, uScale, ...
            Wdes, p, key, Qw, Qu, opts, status);

    otherwise
        error('vtolConstrainedAllocate:UnknownSolver', ...
            'Use ''projectedGN'', ''fmincon'', or ''auto''.');
end

uCommand = uCenter + uScale.*zCommand;
uCommand = min(max(uCommand, lowerBound), upperBound);
[Wactual, wrenchDetail] = vtolWrench(uCommand, p, key);
residual = Wdes-Wactual;

% Utilization: zero at the center of the permissible physical range and one
% at either position limit.  This deliberately captures lower thrust limits,
% unlike T/Tmax alone.  upperFraction is also returned for motor reporting.
utilization = abs((uCommand-uCenter)./uScale);
upperFraction = zeros(size(uCommand));
for i = 1:numel(uCommand)
    if physicalLower(i) >= 0
        upperFraction(i) = uCommand(i)/physicalUpper(i);
    else
        upperFraction(i) = abs(uCommand(i))/max(abs(physicalLower(i)), ...
            abs(physicalUpper(i)));
    end
end

active = localActiveConstraints(uCommand, physicalLower, physicalUpper, ...
    lowerBound, upperBound, opts.constraintTolerance, names);
scaledResidual = Sw*residual;

result = struct();
result.u_cmd = uCommand;
result.W_des = Wdes;
result.W_actual = Wactual;
result.residual = residual;
result.scaledResidual = scaledResidual;
result.scaledResidualNorm = norm(scaledResidual);
result.feasible = result.scaledResidualNorm <= opts.feasibilityTolerance;
result.activeConstraints = active;
result.utilization = utilization;
result.u_max = max(utilization);
result.upperFraction = upperFraction;
result.lowerBound = lowerBound;
result.upperBound = upperBound;
result.physicalLowerBound = physicalLower;
result.physicalUpperBound = physicalUpper;
result.optionsUsed = opts;
result.status = status;
result.wrenchDetail = wrenchDetail;
end

function opts = localOptions(p, userOptions, lowerBound, upperBound)
n = numel(lowerBound);
opts = struct();
opts.solver = localField(userOptions,'solver','projectedGN');
opts.lambdaU = localField(userOptions,'lambdaU',p.allocation.lambdaU);
opts.maxIterations = localField(userOptions,'maxIterations',p.allocation.maxIterations);
opts.dampingInitial = localField(userOptions,'dampingInitial',p.allocation.dampingInitial);
opts.stepTolerance = localField(userOptions,'stepTolerance',p.allocation.stepTolerance);
opts.optimalityTolerance = localField(userOptions,'optimalityTolerance',p.allocation.optimalityTolerance);
opts.constraintTolerance = localField(userOptions,'constraintTolerance',p.allocation.constraintTolerance);
opts.feasibilityTolerance = localField(userOptions,'feasibilityTolerance',p.allocation.feasibilityTolerance);
opts.enforceRateBounds = localField(userOptions,'enforceRateBounds',false);
opts.dt = localField(userOptions,'dt',NaN);
opts.fminconAlgorithm = localField(userOptions,'fminconAlgorithm','sqp');
opts.backtrackFactor = localField(userOptions,'backtrackFactor',p.allocation.backtrackFactor);
opts.minimumStep = localField(userOptions,'minimumStep',p.allocation.minimumStep);
opts.dampingIncrease = localField(userOptions,'dampingIncrease',p.allocation.dampingIncrease);
opts.dampingDecrease = localField(userOptions,'dampingDecrease',p.allocation.dampingDecrease);
opts.maximumDamping = localField(userOptions,'maximumDamping',p.allocation.maximumDamping);

defaultSw = diag(p.allocation.priority(:)./p.allocation.wrenchScale(:));
if isfield(userOptions,'Sw') && ~isempty(userOptions.Sw)
    opts.Sw = userOptions.Sw;
else
    opts.Sw = defaultSw;
end
if isfield(userOptions,'Su') && ~isempty(userOptions.Su)
    opts.Su = userOptions.Su;
else
    opts.Su = eye(n);
end

if ~isequal(size(opts.Sw),[6 6]) || ~isequal(size(opts.Su),[n n])
    error('vtolConstrainedAllocate:WeightSize', ...
        'Sw must be 6x6 and Su must match the actuator vector.');
end
end

function [lowerBound, upperBound] = localEffectiveBounds( ...
    physicalLower, physicalUpper, rateBound, uPrev, opts)
lowerBound = physicalLower;
upperBound = physicalUpper;
if opts.enforceRateBounds && isfinite(opts.dt) && opts.dt > 0
    lowerBound = max(lowerBound, uPrev-rateBound*opts.dt);
    upperBound = min(upperBound, uPrev+rateBound*opts.dt);
end
end

function [z, status] = localProjectedGaussNewton( ...
    z, zLower, zUpper, zPrev, uCenter, uScale, Wdes, p, key, Qw, Qu, opts, status)

z = min(max(z, zLower), zUpper);
damping = opts.dampingInitial;
status.solver = 'projectedGN';

for k = 1:opts.maxIterations
    [objective, gradient, hessian] = localObjective( ...
        z, zPrev, uCenter, uScale, Wdes, p, key, Qw, Qu, opts.lambdaU);
    projectedGradient = z - min(max(z-gradient, zLower), zUpper);
    if norm(projectedGradient, inf) <= opts.optimalityTolerance
        status.exitflag = 1;
        status.iterations = k;
        status.message = 'Projected-gradient tolerance met.';
        return;
    end

    step = -(hessian + damping*eye(numel(z)))\gradient;
    if norm(step, inf) <= opts.stepTolerance
        status.exitflag = 2;
        status.iterations = k;
        status.message = 'Gauss-Newton step tolerance met.';
        return;
    end

    accepted = false;
    lineScale = 1.0;
    while lineScale >= opts.minimumStep
        candidate = min(max(z + lineScale*step, zLower), zUpper);
        if norm(candidate-z, inf) <= opts.stepTolerance
            lineScale = lineScale*opts.backtrackFactor;
            continue;
        end
        candidateObjective = localObjective( ...
            candidate, zPrev, uCenter, uScale, Wdes, p, key, Qw, Qu, opts.lambdaU);
        if candidateObjective < objective
            z = candidate;
            damping = max(damping*opts.dampingDecrease, eps);
            accepted = true;
            break;
        end
        lineScale = lineScale*opts.backtrackFactor;
    end

    if ~accepted
        damping = damping*opts.dampingIncrease;
        if damping > opts.maximumDamping
            status.exitflag = -2;
            status.iterations = k;
            status.message = 'Damping limit reached before convergence.';
            return;
        end
    end
end

status.exitflag = 0;
status.iterations = opts.maxIterations;
status.message = 'Maximum iterations reached.';
end

function [z, status] = localFmincon( ...
    zInitial, zLower, zUpper, zPrev, uCenter, uScale, Wdes, p, key, Qw, Qu, opts, status)

objective = @(z)localObjective(z, zPrev, uCenter, uScale, Wdes, ...
    p, key, Qw, Qu, opts.lambdaU);
fopts = optimoptions('fmincon', 'Algorithm', opts.fminconAlgorithm, ...
    'Display', 'none', 'MaxIterations', opts.maxIterations, ...
    'OptimalityTolerance', opts.optimalityTolerance, ...
    'StepTolerance', opts.stepTolerance);
try
    fopts.SpecifyObjectiveGradient = true;
catch
    % Older MATLAB releases simply use finite-difference objective gradients.
end
[z, ~, exitflag, output] = fmincon(objective, zInitial, [], [], [], [], ...
    zLower, zUpper, [], fopts);
status.solver = 'fmincon';
status.exitflag = exitflag;
status.iterations = output.iterations;
status.message = output.message;
end

function [objective, gradient, hessian] = localObjective( ...
    z, zPrev, uCenter, uScale, Wdes, p, key, Qw, Qu, lambdaU)
u = uCenter + uScale.*z;
W = vtolWrench(u, p, key);
J = vtolWrenchJacobian(u, p, key, 'analytic');
Jz = J*diag(uScale);
errorWrench = W-Wdes;
deltaZ = z-zPrev;
objective = 0.5*(errorWrench.'*Qw*errorWrench + ...
    lambdaU*deltaZ.'*Qu*deltaZ);

if nargout > 1
    gradient = Jz.'*Qw*errorWrench + lambdaU*Qu*deltaZ;
end
if nargout > 2
    hessian = Jz.'*Qw*Jz + lambdaU*Qu;
end
end

function active = localActiveConstraints(u, physicalLower, physicalUpper, ...
    lowerBound, upperBound, tolerance, names)
active = struct();
active.names = names;
active.positionLower = abs(u-physicalLower) <= tolerance;
active.positionUpper = abs(u-physicalUpper) <= tolerance;
active.rateLower = abs(u-lowerBound) <= tolerance & lowerBound > physicalLower+tolerance;
active.rateUpper = abs(u-upperBound) <= tolerance & upperBound < physicalUpper-tolerance;
active.any = active.positionLower | active.positionUpper | ...
    active.rateLower | active.rateUpper;
end

function value = localField(s, name, defaultValue)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end
