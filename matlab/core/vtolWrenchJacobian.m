function J = vtolWrenchJacobian(u, p, architecture, method)
%VTOLWRENCHJACOBIAN Exact or finite-difference dW/du for vtolWrench.
%
%   J = vtolWrenchJacobian(u,p,architecture) uses the analytic derivative.
%   J = vtolWrenchJacobian(...,'central') uses a bound-aware central finite
%   difference for regression testing.  Thrust perturbations and angular
%   perturbations come from p.numerics, not from a hidden constant.

if nargin < 4 || isempty(method)
    method = 'analytic';
end

key = vtolArchitectureKey(architecture);
[lowerBound, upperBound, ~, ~] = vtolActuatorBounds(p, key);
u = min(max(u(:), lowerBound), upperBound);

switch lower(char(method))
    case {'analytic','exact'}
        J = localAnalyticJacobian(u, p, key);
    case {'central','finite-difference','fd'}
        J = localFiniteDifferenceJacobian(u, p, key, lowerBound, upperBound);
    otherwise
        error('vtolWrenchJacobian:UnknownMethod', ...
            'Use ''analytic'' or ''central''.');
end
end

function J = localAnalyticJacobian(u, p, key)
nActuator = numel(u);
J = zeros(6, nActuator);
mainIndices = [1 2 3; 4 5 6];

for i = 1:2
    index = mainIndices(i,:);
    thrust = u(index(1));
    alpha = u(index(2));
    beta = u(index(3));
    ca = cos(alpha); sa = sin(alpha);
    cb = cos(beta);  sb = sin(beta);

    n = [sa; ca*sb; -ca*cb];
    dnDalpha = [ca; -sa*sb; sa*cb];
    dnDbeta = [0.0; ca*cb; ca*sb];
    dF = [n, thrust*dnDalpha, thrust*dnDbeta];
    A = localSkew(p.rMain(:,i)) - p.sigmaMain(i)*p.kappa*eye(3);
    J(:,index) = [dF; A*dF];
end

if strcmp(key, 'tail')
    thrust = u(7);
    delta = u(8);
    n = [-sin(delta); 0.0; -cos(delta)];
    dnDdelta = [-cos(delta); 0.0; sin(delta)];
    dF = [n, thrust*dnDdelta];
    A = localSkew(p.rTail) - p.sigmaTail*p.kappa*eye(3);
    J(:,7:8) = [dF; A*dF];
end
end

function J = localFiniteDifferenceJacobian(u, p, key, lowerBound, upperBound)
nActuator = numel(u);
J = zeros(6, nActuator);
for j = 1:nActuator
    if mod(j-1, 3) == 0 && j <= 6
        h = p.numerics.fdStepThrust;
    elseif strcmp(key,'tail') && j == 7
        h = p.numerics.fdStepThrust;
    else
        h = p.numerics.fdStepAngle;
    end

    uPlus = u;
    uMinus = u;
    uPlus(j) = min(u(j)+h, upperBound(j));
    uMinus(j) = max(u(j)-h, lowerBound(j));
    denominator = uPlus(j)-uMinus(j);

    if denominator <= eps
        J(:,j) = 0.0;
    else
        J(:,j) = (vtolWrench(uPlus,p,key) - vtolWrench(uMinus,p,key)) / denominator;
    end
end
end

function S = localSkew(r)
S = [0.0, -r(3), r(2); ...
     r(3), 0.0, -r(1); ...
    -r(2), r(1), 0.0];
end
