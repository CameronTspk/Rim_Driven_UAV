function [uTrim, Wtrim, info] = vtolHoverTrim(p, architecture)
%VTOLHOVERTRIM Find an exact bounded zero-moment hover trim.
%
%   [uTrim,Wtrim] = vtolHoverTrim(p,architecture) solves
%   W = [0 0 -m*g 0 0 0].' with the same nonlinear model and actuator
%   bounds used in allocation.  The tail configuration includes a small
%   differential-main-alpha seed to cancel the tail reaction torque; the
%   final allocation removes all remaining nonlinear trim residual.

key = vtolArchitectureKey(architecture);
[lowerBound, upperBound, ~, ~] = vtolActuatorBounds(p, key);
Wdesired = [0; 0; -p.m*p.g; 0; 0; 0];

switch key
    case 'twin'
        mainThrust = 0.5*p.m*p.g;
        uSeed = [mainThrust; 0; 0; mainThrust; 0; 0];

    case 'tail'
        mainX = mean(p.rMain(1,:));
        tailArm = abs(p.rTail(1));
        if mainX > 0 && tailArm > 0
            mainThrust = p.m*p.g/(2 + 2*mainX/tailArm);
            tailThrust = 2*mainX*mainThrust/tailArm;
        else
            tailThrust = p.tail.trimThrustGuess;
            mainThrust = 0.5*(p.m*p.g-tailThrust);
        end
        tailThrust = min(max(tailThrust,p.tail.Tmin),p.tail.Tmax);
        mainThrust = min(max(mainThrust,p.main.Tmin),p.main.Tmax);

        % At delta=0, the tail yaw reaction torque is +sigma*kappa*T3.
        % alpha1=+alpha, alpha2=-alpha generates approximately
        % -2*d*Tmain*sin(alpha) in Mz, with no first-order net Fx.
        lateralArm = abs(p.rMain(2,1));
        alphaSeed = 0.0;
        denominator = 2*lateralArm*mainThrust;
        if denominator > 0
            argument = p.sigmaTail*p.kappa*tailThrust/denominator;
            alphaSeed = asin(min(max(argument,-1.0),1.0));
        end
        uSeed = [mainThrust; alphaSeed; 0; ...
                 mainThrust; -alphaSeed; 0; ...
                 tailThrust; 0];
end

uSeed = min(max(uSeed,lowerBound),upperBound);
trimOptions = struct();
trimOptions.solver = 'projectedGN';
trimOptions.Sw = diag(1./p.scaling.wrench);
trimOptions.lambdaU = p.trim.lambdaU;
trimOptions.maxIterations = p.trim.maxIterations;
trimOptions.feasibilityTolerance = p.trim.feasibilityTolerance;
trimOptions.enforceRateBounds = false;
trimOptions.uInitial = uSeed;

allocation = vtolConstrainedAllocate(Wdesired, uSeed, p, key, trimOptions);
uTrim = allocation.u_cmd;
Wtrim = allocation.W_actual;

info = allocation;
info.seed = uSeed;
info.trimResidual = Wdesired-Wtrim;
info.isTrimmed = norm(info.trimResidual./p.scaling.wrench) <= ...
    p.trim.feasibilityTolerance;
end
