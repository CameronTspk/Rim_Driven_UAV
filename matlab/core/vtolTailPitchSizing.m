function sizing = vtolTailPitchSizing(MyIncrementRequired, p, tailTrimThrust)
%VTOLTAILPITCHSIZING Static tail-only pitch-authority screening calculation.
%
%   sizing = vtolTailPitchSizing(MyIncrementRequired,p,T3trim) reports the
%   tail thrust stroke and arm needed to generate +/-MyIncrementRequired at
%   delta=0.  It is a necessary, conservative screening calculation, not a
%   substitute for the 6-DOF combined-maneuver gate because the main props
%   and actuator dynamics also contribute in the full plant.

if ~isfield(p,'tail') || ~isfield(p.tail,'enabled') || ~p.tail.enabled
    error('vtolTailPitchSizing:NoTail', ...
        'Tail-only pitch sizing is defined only for the tail architecture.');
end
if nargin < 3 || isempty(tailTrimThrust)
    tailTrimThrust = p.tail.trimThrustGuess;
end

tailArm = abs(p.rTail(1));
requiredMagnitude = abs(MyIncrementRequired);
increaseMargin = p.tail.Tmax-tailTrimThrust;
decreaseMargin = tailTrimThrust-p.tail.Tmin;

% dMy/dT3 = -tailArm at delta=0.  Hence a negative pitch increment uses
% an increase in tail thrust, while a positive increment uses a decrease.
sizing = struct();
sizing.MyIncrementRequired = requiredMagnitude;
sizing.tailArm = tailArm;
sizing.tailTrimThrust = tailTrimThrust;
sizing.deltaMyNegativeAvailable = tailArm*increaseMargin;
sizing.deltaMyPositiveAvailable = tailArm*decreaseMargin;
sizing.deltaTRequiredEachSign = requiredMagnitude/tailArm;
sizing.minimumArmNegativeMy = requiredMagnitude/max(increaseMargin,eps);
sizing.minimumArmPositiveMy = requiredMagnitude/max(decreaseMargin,eps);
sizing.minimumArmBothSigns = requiredMagnitude/max(min(increaseMargin,decreaseMargin),eps);
sizing.requiredTailThrustStrokeBothSigns = 2*sizing.deltaTRequiredEachSign;
sizing.currentRangeSupportsBothSigns = ...
    (increaseMargin >= sizing.deltaTRequiredEachSign) && ...
    (decreaseMargin >= sizing.deltaTRequiredEachSign);
end
