function [wrench, detail] = vtolWrench(u, p, architecture)
%VTOLWRENCH Exact nonlinear actuator-to-body-wrench mapping.
%
%   wrench = vtolWrench(u,p,architecture) returns
%   [Fx Fy Fz Mx My Mz].' in body axes x-forward, y-right, z-down.
%   The input is projected to the declared physical position limits, so the
%   returned wrench is always physically achievable by position saturation.

key = vtolArchitectureKey(architecture);
[lowerBound, upperBound, ~, names] = vtolActuatorBounds(p, key);
uInput = u(:);

if numel(uInput) ~= numel(lowerBound)
    error('vtolWrench:ActuatorLength', ...
        'Expected %d actuator commands for %s architecture.', ...
        numel(lowerBound), key);
end

uApplied = min(max(uInput, lowerBound), upperBound);
positionSaturated = abs(uApplied-uInput) > p.numerics.boundTolerance;

force = zeros(3,1);
moment = zeros(3,1);
propulsor = repmat(struct('r',zeros(3,1),'n',zeros(3,1), ...
    'F',zeros(3,1),'Q',zeros(3,1),'M',zeros(3,1)), 1, 2 + strcmp(key,'tail'));

% Main propulsors: n = [sin(alpha); cos(alpha)sin(beta); -cos(alpha)cos(beta)].
mainIndices = [1 2 3; 4 5 6];
for i = 1:2
    index = mainIndices(i,:);
    thrust = uApplied(index(1));
    alpha = uApplied(index(2));
    beta = uApplied(index(3));
    direction = [sin(alpha); ...
                 cos(alpha)*sin(beta); ...
                -cos(alpha)*cos(beta)];
    Fi = thrust*direction;
    Qi = -p.sigmaMain(i)*p.kappa*Fi;
    Mi = cross(p.rMain(:,i), Fi) + Qi;

    force = force + Fi;
    moment = moment + Mi;
    propulsor(i) = struct('r',p.rMain(:,i),'n',direction, ...
        'F',Fi,'Q',Qi,'M',Mi);
end

if strcmp(key, 'tail')
    tailThrust = uApplied(7);
    delta = uApplied(8);
    direction = [-sin(delta); 0.0; -cos(delta)];
    Fi = tailThrust*direction;
    Qi = -p.sigmaTail*p.kappa*Fi;
    Mi = cross(p.rTail, Fi) + Qi;

    force = force + Fi;
    moment = moment + Mi;
    propulsor(3) = struct('r',p.rTail,'n',direction, ...
        'F',Fi,'Q',Qi,'M',Mi);
end

wrench = [force; moment];

if nargout > 1
    detail = struct();
    detail.uInput = uInput;
    detail.uApplied = uApplied;
    detail.positionSaturated = positionSaturated;
    detail.actuatorNames = names;
    detail.propulsor = propulsor;
    detail.force = force;
    detail.moment = moment;
end
end
