function metrics = vtolScaleJacobian(J, p, architecture)
%VTOLSCALEJACOBIAN Dimensionless SVD/rank metrics for an effectiveness map.
%
%   Jscaled = diag(1./Wscale) * J * diag(Uscale), where Uscale is half
%   actuator travel and Wscale is p.scaling.wrench.  This avoids assigning
%   physical meaning to an unscaled condition number containing N and rad.

[lowerBound, upperBound, ~, names] = vtolActuatorBounds(p, architecture);
actuatorScale = 0.5*(upperBound-lowerBound);
wrenchScale = p.scaling.wrench(:);

Jscaled = diag(1./wrenchScale) * J * diag(actuatorScale);
[~, S, ~] = svd(Jscaled, 'econ');
singularValues = diag(S);

if isempty(singularValues)
    rankValue = 0;
    conditionNumber = Inf;
    threshold = NaN;
else
    threshold = p.numerics.rankRelativeTolerance*singularValues(1);
    rankValue = sum(singularValues > threshold);
    if rankValue < min(size(Jscaled)) || singularValues(end) <= threshold
        conditionNumber = Inf;
    else
        conditionNumber = singularValues(1)/singularValues(end);
    end
end

metrics = struct();
metrics.J = J;
metrics.Jscaled = Jscaled;
metrics.actuatorScale = actuatorScale;
metrics.wrenchScale = wrenchScale;
metrics.actuatorNames = names;
metrics.singularValues = singularValues;
metrics.rank = rankValue;
metrics.rankThreshold = threshold;
metrics.conditionNumber = conditionNumber;
end
