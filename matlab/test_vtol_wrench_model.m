function report = test_vtol_wrench_model()
%TEST_VTOL_WRENCH_MODEL Regression tests for the exact model and allocator.
%
%   report = test_vtol_wrench_model() runs without a Simulink model.  It
%   throws an error if any test fails and returns all scalar error values for
%   logging in a CI or parameter-sweep workflow.

pTwin = vtolDefaultParams('twin');
pTail = vtolDefaultParams('tail');
[uTwin, WTwin] = vtolHoverTrim(pTwin,'twin');
[uTail, WTail] = vtolHoverTrim(pTail,'tail');

checks = struct('name',{},'passed',{},'value',{},'tolerance',{});

% 1. Exact hover trims.
targetTwin = [0;0;-pTwin.m*pTwin.g;0;0;0];
targetTail = [0;0;-pTail.m*pTail.g;0;0;0];
checks(end+1) = localCheck('Twin hover trim', ...
    norm(WTwin-targetTwin,inf), 1e-5);
checks(end+1) = localCheck('Tail hover trim', ...
    norm(WTail-targetTail,inf), 1e-5);

% 2. Symmetric main alpha generates the expected positive Fx.
alphaTest = deg2rad(10.0);
u = uTwin;
u([2 5]) = alphaTest;
W = vtolWrench(u,pTwin,'twin');
expectedFx = (u(1)+u(4))*sin(alphaTest);
checks(end+1) = localCheck('Symmetric alpha -> Fx', ...
    abs(W(1)-expectedFx), 1e-10);

% 3. Differential beta produces the expected signed weak pitch moment.
betaTest = deg2rad(10.0);
u = uTwin;
u(3) = betaTest;
u(6) = -betaTest;
W = vtolWrench(u,pTwin,'twin');
expectedMy = -pTwin.kappa*(u(1)+u(4))*sin(betaTest);
checks(end+1) = localCheck('Differential beta -> signed My', ...
    abs(W(5)-expectedMy), 1e-10);

% 4. Tail thrust has dMy/dT3=-l_t at delta=0.
u = uTail;
u(8) = 0.0;
thrustStep = 1e-4;
uPlus = u;  uPlus(7) = min(u(7)+thrustStep,pTail.tail.Tmax);
uMinus = u; uMinus(7) = max(u(7)-thrustStep,pTail.tail.Tmin);
dMyDtailThrust = (vtolWrench(uPlus,pTail,'tail') - ...
    vtolWrench(uMinus,pTail,'tail'))/(uPlus(7)-uMinus(7));
checks(end+1) = localCheck('Tail thrust -> signed My', ...
    abs(dMyDtailThrust(5)+abs(pTail.rTail(1))), 1e-7);

% 5. With r3z=0, tail tilt has zero first-order pitch effect at delta=0.
u = uTail;
u(8) = 0.0;
deltaStep = pTail.numerics.fdStepAngle;
uPlus = u;  uPlus(8) = deltaStep;
uMinus = u; uMinus(8) = -deltaStep;
dMyDdelta = (vtolWrench(uPlus,pTail,'tail') - ...
    vtolWrench(uMinus,pTail,'tail'))/(2*deltaStep);
checks(end+1) = localCheck('Tail tilt first-order My at h=0', ...
    abs(dMyDdelta(5)), 1e-7);

% 6. Analytical and bound-aware central finite-difference Jacobians agree.
for item = 1:2
    if item == 1
        p = pTwin; key = 'twin'; u = uTwin;
    else
        p = pTail; key = 'tail'; u = uTail;
    end
    Janalytic = vtolWrenchJacobian(u,p,key,'analytic');
    Jfinite = vtolWrenchJacobian(u,p,key,'central');
    checks(end+1) = localCheck(['Jacobian agreement: ' key], ...
        max(abs(Janalytic(:)-Jfinite(:))), 1e-5);
end

% 7. Position and rate constraints are both obeyed.
[lowerBound, upperBound, rateBound] = vtolActuatorBounds(pTwin,'twin');
rateOptions = struct('solver','projectedGN','enforceRateBounds',true, ...
    'dt',0.002,'maxIterations',40);
hardWrench = [100; 0; -pTwin.m*pTwin.g; 0; 10; 0];
rateResult = vtolConstrainedAllocate(hardWrench,uTwin,pTwin,'twin',rateOptions);
boundViolation = max([lowerBound-rateResult.u_cmd; ...
    rateResult.u_cmd-upperBound; 0]);
rateViolation = max(abs(rateResult.u_cmd-uTwin)-rateBound*rateOptions.dt);
checks(end+1) = localCheck('Allocator position bounds', boundViolation, 1e-10);
checks(end+1) = localCheck('Allocator rate bounds', rateViolation, 1e-10);

% 8. An infeasible wrench returns a residual rather than exceeding a limit.
infeasibleOptions = struct('solver','projectedGN','maxIterations',100);
infeasibleResult = vtolConstrainedAllocate(hardWrench,uTwin,pTwin, ...
    'twin',infeasibleOptions);
checks(end+1) = localCheck('Infeasible wrench leaves residual', ...
    double(infeasibleResult.feasible), 0.0);

report = struct();
report.checks = checks;
report.allPassed = all([checks.passed]);
report.uTwinTrim = uTwin;
report.uTailTrim = uTail;
report.WTwinTrim = WTwin;
report.WTailTrim = WTail;
report.rateLimitedResult = rateResult;
report.infeasibleResult = infeasibleResult;

if nargout == 0 || ~report.allPassed
    for k = 1:numel(checks)
        fprintf('%-42s %s  value=% .3e  tol=% .3e\n', checks(k).name, ...
            localPassFail(checks(k).passed), checks(k).value, checks(k).tolerance);
    end
end
if ~report.allPassed
    error('test_vtol_wrench_model:Failure', ...
        'One or more nonlinear-wrench regression tests failed.');
end
end

function check = localCheck(name, value, tolerance)
check = struct();
check.name = name;
check.value = value;
check.tolerance = tolerance;
check.passed = value <= tolerance;
end

function label = localPassFail(passed)
if passed
    label = 'PASS';
else
    label = 'FAIL';
end
end
