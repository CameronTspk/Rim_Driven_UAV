function summary = plotCombinedWrenchDiagnostic(logData, p, architecture, outputDirectory)
%PLOTCOMBINEDWRENCHDIAGNOSTIC Plot allocation and actuator losses in 6-DOF logs.
%
% Required numeric logData fields (one row per time step):
%   t, W_cmd, W_actual, u_actual
% Optional fields:
%   W_allocated  -- wrench immediately after static allocation, before lag
%   u_cmd        -- static allocator command, before actuator dynamics
%
% This function separates W_cmd-W_allocated (instantaneous constrained
% allocation loss) from W_allocated-W_actual (motor/gimbal rate and lag loss).
% It is deliberately independent of a particular Simulink logging format;
% adapt the signal extraction at the model boundary, then call this function.

if nargin < 4 || isempty(outputDirectory)
    outputDirectory = pwd;
end
if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

key = vtolArchitectureKey(architecture);
[lowerBound, upperBound, ~, names] = vtolActuatorBounds(p,key);
t = logData.t(:);
Wcmd = logData.W_cmd;
Wactual = logData.W_actual;
uActual = logData.u_actual;

localValidateSize(Wcmd,t,6,'W_cmd');
localValidateSize(Wactual,t,6,'W_actual');
localValidateSize(uActual,t,numel(lowerBound),'u_actual');

if isfield(logData,'W_allocated') && ~isempty(logData.W_allocated)
    Wallocated = logData.W_allocated;
    localValidateSize(Wallocated,t,6,'W_allocated');
else
    Wallocated = Wactual;
end
if isfield(logData,'u_cmd') && ~isempty(logData.u_cmd)
    uCommand = logData.u_cmd;
    localValidateSize(uCommand,t,numel(lowerBound),'u_cmd');
else
    uCommand = uActual;
end

uCenter = 0.5*(lowerBound+upperBound).';
uScale = 0.5*(upperBound-lowerBound).';
utilization = abs((uActual-uCenter)./uScale);
uMax = max(utilization,[],2);
allocationLoss = Wcmd-Wallocated;
lagLoss = Wallocated-Wactual;

figure('Color','w','Name','Wrench tracking');
subplot(2,1,1);
plot(t,Wcmd(:,1),'k--','LineWidth',1.1); hold on;
plot(t,Wallocated(:,1),'b','LineWidth',1.1);
plot(t,Wactual(:,1),'r','LineWidth',1.1); grid on;
ylabel('F_x [N]'); legend('commanded','post-allocation','actual','Location','best');
title('Longitudinal-force tracking');
subplot(2,1,2);
plot(t,Wcmd(:,5),'k--','LineWidth',1.1); hold on;
plot(t,Wallocated(:,5),'b','LineWidth',1.1);
plot(t,Wactual(:,5),'r','LineWidth',1.1); grid on;
xlabel('Time [s]'); ylabel('M_y [N m]');
legend('commanded','post-allocation','actual','Location','best');
title('Pitch-moment tracking');
saveas(gcf,fullfile(outputDirectory,[key '_wrench_tracking.png']));

figure('Color','w','Name','Main gimbals');
gimbalIndex = [2 3 5 6];
for k = 1:numel(gimbalIndex)
    subplot(2,2,k);
    index = gimbalIndex(k);
    plot(t,rad2deg(uCommand(:,index)),'--','LineWidth',1.0); hold on;
    plot(t,rad2deg(uActual(:,index)),'LineWidth',1.1);
    yline(rad2deg(upperBound(index)),'r:');
    yline(rad2deg(lowerBound(index)),'r:'); grid on;
    xlabel('Time [s]'); ylabel([names{index} ' [deg]']);
    title(names{index});
end
legend('command','actual','limit','Location','best');
saveas(gcf,fullfile(outputDirectory,[key '_main_gimbals.png']));

figure('Color','w','Name','Thrust and utilization');
subplot(2,1,1);
thrustIndex = [1 4];
if strcmp(key,'tail')
    thrustIndex = [thrustIndex 7];
end
plot(t,uActual(:,thrustIndex),'LineWidth',1.1); hold on;
for k = 1:numel(thrustIndex)
    yline(upperBound(thrustIndex(k)),'k:');
end
grid on; ylabel('Thrust [N]');
legend(names(thrustIndex),'Location','best'); title('Individual propeller thrusts');
subplot(2,1,2);
plot(t,uMax,'k','LineWidth',1.2); hold on; yline(0.9,'r:'); yline(1.0,'r--');
grid on; xlabel('Time [s]'); ylabel('u_{max}');
title('Normalized position-limit utilization (one = a position limit)');
saveas(gcf,fullfile(outputDirectory,[key '_thrust_and_utilization.png']));

summary = struct();
summary.peakFxAllocationLoss = max(abs(allocationLoss(:,1)));
summary.peakMyAllocationLoss = max(abs(allocationLoss(:,5)));
summary.peakFxLagLoss = max(abs(lagLoss(:,1)));
summary.peakMyLagLoss = max(abs(lagLoss(:,5)));
summary.peakUtilization = max(uMax);
summary.saturationTime = trapz(t,double(uMax >= 0.999));
summary.highUtilizationFraction = mean(uMax >= 0.9);
summary.utilization = utilization;
summary.uMax = uMax;
summary.allocationLoss = allocationLoss;
summary.lagLoss = lagLoss;
end

function localValidateSize(array,t,nColumn,name)
if size(array,1) ~= numel(t) || size(array,2) ~= nColumn
    error('plotCombinedWrenchDiagnostic:LogSize', ...
        '%s must be N-by-%d, where N=numel(t).',name,nColumn);
end
end
