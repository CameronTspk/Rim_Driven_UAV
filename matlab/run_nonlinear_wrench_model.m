%% RUN NONLINEAR VTOL WRENCH MODEL — NO FOLDER PICKER

close all;
clc;

% ===== CHANGE ONLY THIS LINE =====
modelDir = 'C:\Users\cthom\OneDrive\Desktop\Classes\Git_substance\nonlinear_wrench_model_final\nonlinear_wrench_model';
% =================================

runFullSweep = false;          % false = quick run; true = denser/slower run
peakPitchMomentScreen = 2.682; % N m; prior combined-test screening value

% Set MATLAB's working folder and add all model functions to its path.
assert(exist(modelDir,'dir') == 7, ...
    'Folder does not exist. Check modelDir.');

cd(modelDir);
addpath(modelDir);

% Confirm required files exist.
requiredFiles = {'vtolDefaultParams.m','vtolWrench.m', ...
    'vtolConstrainedAllocate.m','vtolHoverTrim.m', ...
    'vtolFeasibleWrenchMap.m','test_vtol_wrench_model.m', ...
    'run_vtol_wrench_analysis.m'};

for k = 1:numel(requiredFiles)
    assert(exist(requiredFiles{k},'file') == 2, ...
        'Missing %s. Check that modelDir is the nonlinear_wrench_model folder.', ...
        requiredFiles{k});
end

fprintf('\n=== Nonlinear VTOL Wrench Model ===\n');
fprintf('Working folder:\n%s\n\n',modelDir);

% 1. Run regression tests.
fprintf('Running regression tests...\n');

testReport = test_vtol_wrench_model();

assert(testReport.allPassed, ...
    'Regression tests failed. Stop here and send me the MATLAB error.');

fprintf('PASS: all %d regression tests passed.\n\n', ...
    numel(testReport.checks));

% 2. Generate constrained nonlinear authority maps and plots.
fprintf('Running nonlinear authority-map analysis...\n');

results = run_vtol_wrench_analysis( ...
    runFullSweep,peakPitchMomentScreen);

% 3. Demonstrate a single allocation request for each architecture.
pTwin = vtolDefaultParams('twin');
[uTwinTrim,~] = vtolHoverTrim(pTwin,'twin');

WtwinRequest = [2.0; 0.0; -pTwin.m*pTwin.g; 0.0; 0.10; 0.0];

twinAllocation = vtolConstrainedAllocate( ...
    WtwinRequest,uTwinTrim,pTwin,'twin', ...
    struct('solver','projectedGN'));

pTail = vtolDefaultParams('tail');
[uTailTrim,~] = vtolHoverTrim(pTail,'tail');

WtailRequest = [2.0; 0.0; -pTail.m*pTail.g; 0.0; 0.10; 0.0];

tailAllocation = vtolConstrainedAllocate( ...
    WtailRequest,uTailTrim,pTail,'tail', ...
    struct('solver','projectedGN'));

fprintf('\nTwin-only example:\n');
fprintf('  u_max = %.3f\n',twinAllocation.u_max);
fprintf('  scaled residual = %.4g\n',twinAllocation.scaledResidualNorm);

fprintf('\nTail-assisted example:\n');
fprintf('  u_max = %.3f\n',tailAllocation.u_max);
fprintf('  scaled residual = %.4g\n',tailAllocation.scaledResidualNorm);

% Keep outputs available for inspection.
assignin('base','nonlinearWrenchResults',results);
assignin('base','twinAllocationDemo',twinAllocation);
assignin('base','tailAllocationDemo',tailAllocation);

fprintf('\nFinished.\n');
fprintf('Open results here:\n%s\n',fullfile(modelDir,'results'));