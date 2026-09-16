function p = vtolDefaultParams(architecture)
%VTOLDEFAULTPARAMS Reference parameters for the nonlinear wrench model.
%
%   p = vtolDefaultParams('twin') returns the twin-only baseline.
%   p = vtolDefaultParams('tail') returns the tail-assisted prior-test
%   baseline.  All values are SI.  Replace geometry, bench-test, and CAD
%   values here rather than changing the model equations.

key = vtolArchitectureKey(architecture);

p = struct();
p.architecture = key;
p.g = 9.81;
p.rho = 1.225;

% Rim-driven propulsor placeholder coefficients.  kappa is the reaction
% torque arm in Q = kappa*T, with spin sign carried separately.
p.prop.D = 0.180;
p.prop.CT = 0.10;
p.prop.CQ = 0.008;
p.kappa = (p.prop.CQ / p.prop.CT) * p.prop.D;

% Main-prop actuator and rate limits.
p.main.Tmin = 1.0;
p.main.Tmax = 14.0;
p.main.alphaMax = deg2rad(35.0);
p.main.betaMax = deg2rad(35.0);
p.main.TdotMax = 200.0;
p.main.alphaDotMax = deg2rad(300.0);
p.main.betaDotMax = deg2rad(300.0);
p.main.tauMotor = 0.030;
p.main.tauGimbal = 0.030;
p.sigmaMain = [1; -1];

switch key
    case 'twin'
        p.m = 1.50;
        p.I = diag([0.070, 0.020, 0.075]);
        p.geometry.mainHalfSpan = 0.300;
        p.rMain = [0.0, 0.0; ...
                   p.geometry.mainHalfSpan, -p.geometry.mainHalfSpan; ...
                   0.0, 0.0];
        p.tail.enabled = false;
        p.rTail = [0; 0; 0];
        p.sigmaTail = 0;

    case 'tail'
        p.m = 1.75;
        p.I = diag([0.070, 0.05514678, 0.11014678]);
        p.geometry.mainHalfSpan = 0.300;
        p.geometry.tailArm = 0.350;
        p.geometry.mainForwardOffset = 0.1063;
        p.rMain = [p.geometry.mainForwardOffset, p.geometry.mainForwardOffset; ...
                   p.geometry.mainHalfSpan, -p.geometry.mainHalfSpan; ...
                   0.0, 0.0];
        p.rTail = [-p.geometry.tailArm; 0.0; 0.0];
        p.tail.enabled = true;
        p.tail.Tmin = 0.5;
        p.tail.Tmax = 10.0;
        p.tail.deltaMax = deg2rad(45.0);
        p.tail.TdotMax = 200.0;
        p.tail.deltaDotMax = deg2rad(300.0);
        p.tail.tauMotor = 0.030;
        p.tail.tauGimbal = 0.030;
        p.tail.trimThrustGuess = 4.0;
        p.sigmaTail = 1;
end

% Numerics used by the exact model and finite-difference regression check.
p.numerics.fdStepThrust = 1e-5;       % N
p.numerics.fdStepAngle = 1e-6;        % rad
p.numerics.rankRelativeTolerance = 1e-8;
p.numerics.boundTolerance = 1e-9;

% Scales used only for dimensionless conditioning.  They are deliberately
% independent of the allocation-priority weights below.
referenceLength = max([sqrt(sum(p.rMain.^2,1)), norm(p.rTail), 1e-3]);
p.scaling.referenceLength = referenceLength;
p.scaling.wrench = [p.m*p.g; p.m*p.g; p.m*p.g; ...
                    p.m*p.g*referenceLength; ...
                    p.m*p.g*referenceLength; ...
                    p.m*p.g*referenceLength];

% Allocation scales and priorities.  Moment priorities exceed translational
% priorities when the requested wrench is infeasible, with extra emphasis on
% the known weak pitch axis.  Tune these only after changing p.scaling.
p.allocation.wrenchScale = [5.0; 5.0; p.m*p.g; 0.50; 0.25; 0.50];
p.allocation.priority = [1.0; 1.0; 1.0; 4.0; 8.0; 4.0];
p.allocation.lambdaU = 1e-4;
p.allocation.maxIterations = 80;
p.allocation.dampingInitial = 1e-4;
p.allocation.stepTolerance = 1e-8;
p.allocation.optimalityTolerance = 1e-7;
p.allocation.constraintTolerance = 1e-7;
p.allocation.feasibilityTolerance = 2e-2;
p.allocation.backtrackFactor = 0.5;
p.allocation.minimumStep = 1e-5;
p.allocation.dampingIncrease = 10.0;
p.allocation.dampingDecrease = 0.3;
p.allocation.maximumDamping = 1e12;

% Exact-trim solve settings.  This corrects the tail-prop reaction-torque
% bias that is present if all gimbal angles are simply set to zero.
p.trim.lambdaU = 1e-8;
p.trim.maxIterations = 200;
p.trim.feasibilityTolerance = 1e-6;

% Map settings.  These are parameter data rather than hidden assumptions in
% the feasibility-map routine.  The driver can override any field.
p.map.nFx = 21;
p.map.nFy = 21;
p.map.nMy = 31;
p.map.nAxis = 41;
p.map.nAngle = 9;
p.map.FxLimit = 2*p.main.Tmax*sin(p.main.alphaMax);
if p.tail.enabled
    p.map.FxLimit = p.map.FxLimit + p.tail.Tmax*sin(p.tail.deltaMax);
end
p.map.FyLimit = 2*p.main.Tmax*sin(p.main.betaMax);
if strcmp(key,'twin')
    % The reaction-torque pitch envelope is O(0.1 N m), so use a narrow,
    % resolved grid rather than burying its boundary in a multi-N m range.
    p.map.nMy = 41;
    p.map.MyLimit = 0.20;
    p.map.MomentLimits = [5.0; 0.20; 5.0];
else
    p.map.MyLimit = 3.5;
    p.map.MomentLimits = [5.0; 3.5; 5.0];
end
p.map.fzTolerance = 0.02*p.m*p.g;
p.map.scaledResidualTolerance = 2e-2;
p.map.solver = 'projectedGN';
p.map.maxIterations = 60;
end
