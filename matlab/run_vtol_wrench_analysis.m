function results = run_vtol_wrench_analysis(fullSweep, peakPitchMomentForScreening)
%RUN_VTOL_WRENCH_ANALYSIS Driver for nonlinear maps, tests, and figures.
%
%   results = run_vtol_wrench_analysis() runs a quick 13-by-13 map suitable
%   for debugging.  run_vtol_wrench_analysis(true) runs the denser default
%   grid defined in p.map.  The second input is a required-pitch-moment
%   screening value for tail-only sizing; the default 2.682 N m is the peak
%   controller request from the prior q(0)=4 rad/s + 1 m x-step study and
%   must be replaced with the logged value after controller integration.

if nargin < 1 || isempty(fullSweep)
    fullSweep = false;
end
if nargin < 2 || isempty(peakPitchMomentForScreening)
    peakPitchMomentForScreening = 2.682;
end

rootDirectory = fileparts(mfilename('fullpath'));
addpath(rootDirectory);
outputDirectory = fullfile(rootDirectory,'results');
if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

fprintf('Running nonlinear wrench regression tests...\n');
testReport = test_vtol_wrench_model();
fprintf('All %d regression tests passed.\n',numel(testReport.checks));

pTwin = vtolDefaultParams('twin');
pTail = vtolDefaultParams('tail');
if fullSweep
    gridSpec = struct();
else
    gridSpec = struct('nFx',13,'nFy',13,'nAxis',25, ...
        'nAngle',7,'maxIterations',50,'solver','projectedGN');
end

fprintf('Computing twin-only constrained wrench map...\n');
mapTwin = vtolFeasibleWrenchMap(pTwin,'twin',gridSpec);
fprintf('Computing tail-assisted constrained wrench map...\n');
mapTail = vtolFeasibleWrenchMap(pTail,'tail',gridSpec);

[uTailTrim,~,~] = vtolHoverTrim(pTail,'tail');
tailSizing = vtolTailPitchSizing(peakPitchMomentForScreening,pTail,uTailTrim(7));

results = struct();
results.testReport = testReport;
results.mapTwin = mapTwin;
results.mapTail = mapTail;
results.tailSizing = tailSizing;
save(fullfile(outputDirectory,'nonlinear_wrench_analysis.mat'),'results','-v7');

localPlotForcePitchEnvelopes(mapTwin,mapTail,outputDirectory);
localPlotGridDiagnostics(mapTwin,mapTail,outputDirectory);
localPlotLateralTrade(mapTwin,mapTail,outputDirectory);
localPlotAuthorityAndConditioning(mapTwin,mapTail,outputDirectory);
localPlotScaledSingularValues(mapTwin,mapTail,outputDirectory);

summary = table( ...
    {'Twin-only';'Tail-assisted'}, ...
    [mapTwin.axisAuthority.angularAccelerationMagnitude(1); mapTail.axisAuthority.angularAccelerationMagnitude(1)], ...
    [mapTwin.axisAuthority.angularAccelerationMagnitude(2); mapTail.axisAuthority.angularAccelerationMagnitude(2)], ...
    [mapTwin.axisAuthority.angularAccelerationMagnitude(3); mapTail.axisAuthority.angularAccelerationMagnitude(3)], ...
    [max(mapTwin.fxMy.uMax(:)); max(mapTail.fxMy.uMax(:))], ...
    [mean(mapTwin.fxMy.feasible(:)); mean(mapTail.fxMy.feasible(:))], ...
    'VariableNames',{'Architecture','MaxRollAccel','MaxPitchAccel','MaxYawAccel', ...
    'PeakUMax','FxMyFeasibleFraction'});
writetable(summary,fullfile(outputDirectory,'nonlinear_wrench_summary.csv'));

fprintf('Results saved in %s\n',outputDirectory);
fprintf(['Tail-only static screening at |delta My| = %.3f N m: current arm %.3f m; ' ...
    'required symmetric tail-only arm %.3f m.\n'], ...
    peakPitchMomentForScreening,tailSizing.tailArm,tailSizing.minimumArmBothSigns);
end

function localPlotForcePitchEnvelopes(mapTwin,mapTail,outputDirectory)
figure('Color','w','Name','Fx and My constrained envelopes');
subplot(2,2,1);
localEnvelopePlot(mapTwin.fxMy,'Twin-only: M_y versus requested F_x','F_x request [N]','M_y achievable [N m]');
subplot(2,2,2);
localEnvelopePlot(mapTail.fxMy,'Tail-assisted: M_y versus requested F_x','F_x request [N]','M_y achievable [N m]');
subplot(2,2,3);
localReverseEnvelopePlot(mapTwin.fxMy,'Twin-only: F_x versus requested M_y','M_y request [N m]','F_x achievable [N]');
subplot(2,2,4);
localReverseEnvelopePlot(mapTail.fxMy,'Tail-assisted: F_x versus requested M_y','M_y request [N m]','F_x achievable [N]');
saveas(gcf,fullfile(outputDirectory,'fx_my_feasible_envelopes.png'));
end

function localPlotGridDiagnostics(mapTwin,mapTail,outputDirectory)
figure('Color','w','Name','Allocation residual and utilization');
subplot(2,2,1);
localGridPlot(mapTwin.fxMy.primaryRequested,mapTwin.fxMy.secondaryRequested, ...
    mapTwin.fxMy.scaledResidualNorm,'Twin-only scaled allocation residual','F_x request [N]','M_y request [N m]');
subplot(2,2,2);
localGridPlot(mapTail.fxMy.primaryRequested,mapTail.fxMy.secondaryRequested, ...
    mapTail.fxMy.scaledResidualNorm,'Tail-assisted scaled allocation residual','F_x request [N]','M_y request [N m]');
subplot(2,2,3);
localGridPlot(mapTwin.fxMy.primaryRequested,mapTwin.fxMy.secondaryRequested, ...
    mapTwin.fxMy.uMax,'Twin-only u_{max}','F_x request [N]','M_y request [N m]');
subplot(2,2,4);
localGridPlot(mapTail.fxMy.primaryRequested,mapTail.fxMy.secondaryRequested, ...
    mapTail.fxMy.uMax,'Tail-assisted u_{max}','F_x request [N]','M_y request [N m]');
saveas(gcf,fullfile(outputDirectory,'allocation_residual_and_utilization.png'));
end

function localPlotLateralTrade(mapTwin,mapTail,outputDirectory)
figure('Color','w','Name','Fy and My trade');
subplot(2,2,1);
localEnvelopePlot(mapTwin.fyMy,'Twin-only: M_y versus requested F_y','F_y request [N]','M_y achievable [N m]');
subplot(2,2,2);
localEnvelopePlot(mapTail.fyMy,'Tail-assisted: M_y versus requested F_y','F_y request [N]','M_y achievable [N m]');
subplot(2,2,3);
localReverseEnvelopePlot(mapTwin.fyMy,'Twin-only: F_y versus requested M_y','M_y request [N m]','F_y achievable [N]');
subplot(2,2,4);
localReverseEnvelopePlot(mapTail.fyMy,'Tail-assisted: F_y versus requested M_y','M_y request [N m]','F_y achievable [N]');
saveas(gcf,fullfile(outputDirectory,'fy_my_feasible_envelopes.png'));
end

function localPlotAuthorityAndConditioning(mapTwin,mapTail,outputDirectory)
figure('Color','w','Name','Angular authority and scaled conditioning');
subplot(2,2,1);
authority = [mapTwin.axisAuthority.angularAccelerationMagnitude, ...
             mapTail.axisAuthority.angularAccelerationMagnitude];
bar(authority); grid on;
set(gca,'XTickLabel',{'M_x','M_y','M_z'});
ylabel('Maximum |angular acceleration| [rad/s^2]');
legend('Twin-only','Tail-assisted','Location','best');
title('Hover-constrained angular authority');

subplot(2,2,2);
localConditionPlot(mapTwin,'Twin-only scaled condition number');
subplot(2,2,3);
localConditionPlot(mapTail,'Tail-assisted scaled condition number');
subplot(2,2,4);
rankData = [mapTwin.jacobian.mode(2).rank(:), mapTail.jacobian.mode(2).rank(:)];
histogram(rankData(:,1),'BinMethod','integers','DisplayStyle','stairs','LineWidth',1.4); hold on;
histogram(rankData(:,2),'BinMethod','integers','DisplayStyle','stairs','LineWidth',1.4);
grid on; xlabel('Scaled Jacobian rank'); ylabel('Grid-count');
legend('Twin-only','Tail-assisted','Location','best');
title('Differential-beta gimbal grid');
saveas(gcf,fullfile(outputDirectory,'authority_and_scaled_conditioning.png'));
end

function localEnvelopePlot(mapData,titleText,xText,yText)
plot(mapData.primaryRequested, ...
    mapData.envelope.secondaryMinimumVsPrimaryRequest, ...
    'b','LineWidth',1.2);
hold on;

plot(mapData.primaryRequested, ...
    mapData.envelope.secondaryMaximumVsPrimaryRequest, ...
    'b','LineWidth',1.2);

yline(0,'k:');
grid on;
xlabel(xText);
ylabel(yText);
title(titleText);
legend('minimum','maximum','Location','best');
end

function localReverseEnvelopePlot(mapData,titleText,xText,yText)
plot(mapData.secondaryRequested, ...
    mapData.envelope.primaryMinimumVsSecondaryRequest, ...
    'b','LineWidth',1.2);
hold on;

plot(mapData.secondaryRequested, ...
    mapData.envelope.primaryMaximumVsSecondaryRequest, ...
    'b','LineWidth',1.2);

yline(0,'k:');
grid on;
xlabel(xText);
ylabel(yText);
title(titleText);
legend('minimum','maximum','Location','best');
end

function localGridPlot(x,y,z,titleText,xText,yText)
imagesc(x,y,z); set(gca,'YDir','normal'); axis tight; colorbar; grid on;
xlabel(xText); ylabel(yText); title(titleText);
end

function localConditionPlot(map,titleText)
condition = map.jacobian.mode(2).conditionNumber;
finiteCondition = condition(isfinite(condition));
if isempty(finiteCondition)
    condition(:) = 1e12;
else
    condition(~isfinite(condition)) = 2*max(finiteCondition);
end
imagesc(rad2deg(map.jacobian.alphaGrid),rad2deg(map.jacobian.betaGrid),log10(condition));
set(gca,'YDir','normal'); axis tight; colorbar; grid on;
xlabel('Common alpha [deg]'); ylabel('Differential beta [deg]');
title(titleText);
end

function localPlotScaledSingularValues(mapTwin,mapTail,outputDirectory)
figure('Color','w','Name','Scaled Jacobian singular values');
subplot(2,2,1);
localMinimumSingularValuePlot(mapTwin,'Twin-only: log_{10}(sigma_{min})');
subplot(2,2,2);
localMinimumSingularValuePlot(mapTail,'Tail-assisted: log_{10}(sigma_{min})');
subplot(2,2,3);
localCenterSpectrumPlot(mapTwin,'Twin-only central spectrum');
subplot(2,2,4);
localCenterSpectrumPlot(mapTail,'Tail-assisted central spectrum');
saveas(gcf,fullfile(outputDirectory,'scaled_jacobian_singular_values.png'));
end

function localMinimumSingularValuePlot(map,titleText)
singular = map.jacobian.mode(2).singularValues;
minimumSingular = squeeze(singular(end,:,:));
imagesc(rad2deg(map.jacobian.alphaGrid),rad2deg(map.jacobian.betaGrid), ...
    log10(max(minimumSingular,eps)));
set(gca,'YDir','normal'); axis tight; colorbar; grid on;
xlabel('Common alpha [deg]'); ylabel('Differential beta [deg]');
title(titleText);
end

function localCenterSpectrumPlot(map,titleText)
ia = ceil(numel(map.jacobian.alphaGrid)/2);
ib = ceil(numel(map.jacobian.betaGrid)/2);
singular = squeeze(map.jacobian.mode(2).singularValues(:,ib,ia));
semilogy(1:numel(singular),singular,'o-','LineWidth',1.2); grid on;
xlabel('Singular-value index'); ylabel('Scaled singular value');
title(titleText);
end
