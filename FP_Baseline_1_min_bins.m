clear all
close all

%% ========================================================================
%% USER CONFIGURATION - MUST FILL OUT
%% ========================================================================

% Directory containing the data
theFolder = '/Users/HL801/Desktop/Photometry_data/CeA_RF_4cup_C4_M/10.26.23-LPSTESTDAY/269437/1L_and_2R/';

% Mouse identifier
mouseName = '2R';

% Sensor type: 'GCaMP' or 'Sensor'
SigProcMode = 'GCaMP';

% Analysis start time (seconds)
beginHere = 0;

% When has the photobleaching period of the baseline ended? 
BASELINE_START_MIN = 0;

% ========================================================================
% DATA FORMAT SELECTION
% ========================================================================
% IMPORTANT: Uncomment ONE of the two sections below based on your data format

% --- OPTION 1: OLD FORMAT (collected before August 15, 2025) ---
dataFormat = 'old';
signalColumn = 'Region3G';      % Options: Region0G, Region1G, Region2G, Region3G
refColumn = 'Region3G';         % Usually same as signal
timeColumn = 'Timestamp';
ledStateColumn = 'LedState';
greenLedState = 6;              % LED state for GCaMP signal; Green is  6 if trig 1, Green is 2 if trig 3
isoLedState = 1;                % LED state for isosbestic

% % --- OPTION 2: NEW FORMAT (collected after August 15, 2025) ---
% dataFormat = 'new';
% signalColumn = 'G1';             % Options: G0, G1, G2, G3
% refColumn = 'G1';                % Usually same as signal
% timeColumn = 'ComputerTimestamp';
% ledStateColumn = 'LedState';
% greenLedState = 2;               % LED state for GCaMP signal
% isoLedState = 1;                 % LED state for isosbestic

% ========================================================================
% EXPERIMENTAL TIMELINE (in seconds)
% ========================================================================

% Injection times
injStart = 1791;
injEnd = 1814;

% Direct social interaction (DSI)
DSIStart = [];
DSIEnd = [];

% Cage addition
CageStart = [];
CageEnd = [];

% Cage with food
CagefoodStart = [];
CagefoodEnd = [];

% Food only
FoodStart = [];
FoodEnd = [];

% Bad periods to remove (artifacts)
rmStart = [];
rmEnd = [];

% Output directory
outFolder = [theFolder '/Final_output_1min-GCaMP-_' mouseName '_' datestr(now, '2025_10_10') '/'];
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

% ========================================================================
% ANALYSIS PARAMETERS
% ========================================================================

% Baseline estimation
BASELINE_WINDOW_MIN = 5;        % Window size for percentile filtering (minutes)
BASELINE_PERCENTILE = 8;        % Percentile for baseline (8th = resistant to transients)

% Transient detection (optimized for GCaMP6s at ~20 Hz)
TRANSIENT_MIN_PROMINENCE = 1.5; % Z-score units                        original input: 1.5
TRANSIENT_MIN_SEPARATION = 2.0; % Seconds (GCaMP6s decay ~0.5-1 sec)   original input: 2.0
TRANSIENT_MIN_WIDTH = 0.5;      % Seconds                              original input: 0.5
TRANSIENT_RISE_MIN = 0.05;      % Minimum rise time (sec)              original input: 0.05
TRANSIENT_RISE_MAX = 20.0;       % Maximum rise time (sec)              original input: 0.8
TRANSIENT_DECAY_MIN = 0.3;      % Minimum decay time (sec)             original input: 0.3
TRANSIENT_DECAY_MAX = 120.0;      % Maximum decay time (sec)             original input: 3.0

% Output downsampling
PROCESSED_SIGNALS_DOWNSAMPLE = 10;  % Factor for downsampling full time series output (1 = no downsampling)

% Binning for output
BIN_SIZE_MIN = 1;               % 1-minute bins for time series output

% Visualization
PLOT_RESULTS = true;

%% ========================================================================
%% LOAD DATA
%% ========================================================================

fprintf('\n========================================\n');
fprintf('LOADING DATA\n');
fprintf('========================================\n');
fprintf('Data format: %s\n', dataFormat);
fprintf('Signal column: %s\n', signalColumn);

cd(theFolder)

% Find photometry file
fpFile = dir([theFolder '/P*.csv']);
if length(fpFile) == 1
    fpFile = fpFile.name;
    fprintf('Found file: %s\n', fpFile);
elseif isempty(fpFile)
    error('No photometry files found in: %s', theFolder);
else
    error('Multiple photometry files found. Please ensure only one CSV file.');
end

% Read the data
fp = readtable(fpFile);

% Display available columns for debugging
fprintf('\nAvailable columns in data file:\n');
disp(fp.Properties.VariableNames);

% Extract signals based on LED state and format
try
    % Get LED states
    ledStates = fp.(ledStateColumn);
    
    % Extract green (GCaMP) signal
    greenMask = (ledStates == greenLedState);
    theSig = fp.(signalColumn)(greenMask);
    theTime = fp.(timeColumn)(greenMask);
    
    % Extract isosbestic reference
    isoMask = (ledStates == isoLedState);
    ref = fp.(refColumn)(isoMask);
    
    % Convert timestamps to seconds (new format uses milliseconds)
    if strcmp(dataFormat, 'new')
        % Check if timestamps are in milliseconds (values > 1000000)
        if mean(theTime) > 1000000
            fprintf('Converting timestamps from milliseconds since epoch to relative seconds...\n');
            % Convert to seconds and zero to start
            theTime = (theTime - theTime(1)) / 1000.0;
            fprintf('Time now spans: 0 to %.2f seconds\n', theTime(end));
        else
            % Just zero the timestamps
            theTime = theTime - theTime(1);
        end
    end

    fprintf('Extracted signals successfully\n');
    fprintf('Green LED frames: %d\n', sum(greenMask));
    fprintf('Isosbestic LED frames: %d\n', sum(isoMask));
    
catch ME
    fprintf('Error extracting data: %s\n', ME.message);
    fprintf('Please check your column names and LED state values\n');
    rethrow(ME);
end

% Ensure equal lengths
theMin = min([length(theSig), length(ref)]);
theSig = theSig(1:theMin);
ref = ref(1:theMin);
theTime = theTime(1:theMin);

% Adjust time to start at zero (if not already done)
if theTime(1) ~= 0
    theTime = theTime - theTime(1);
end

fprintf('Loaded %d frames (%.1f minutes)\n', length(theSig), theTime(end)/60);
fprintf('Effective sampling rate: %.2f Hz\n', length(theSig)/theTime(end));

%% ========================================================================
%% SIGNAL PREPROCESSING
%% ========================================================================

fprintf('\n========================================\n');
fprintf('SIGNAL PREPROCESSING\n');
fprintf('========================================\n');

% Truncate to analysis start point
if beginHere >= theTime(end)
    fprintf('WARNING: beginHere (%.2f) is beyond data end (%.2f). Using full dataset.\n', ...
        beginHere, theTime(end));
    beginHere = 0;
end

beginHereIdx = find(theTime >= beginHere, 1);
if isempty(beginHereIdx)
    beginHereIdx = 1;
end

theSig = theSig(beginHereIdx:end);
ref = ref(beginHereIdx:end);
theTime = theTime(beginHereIdx:end) - theTime(beginHereIdx);

fprintf('Starting analysis at %.2f seconds\n', beginHere);
fprintf('Frames after truncation: %d (%.1f minutes)\n', length(theSig), theTime(end)/60);

% Isosbestic correction (if GCaMP mode)
if strcmp(SigProcMode, 'GCaMP')
    fprintf('Performing isosbestic correction...\n');
    
    % Fit exponential to reference signal (bleaching correction)
    temp_x = (1:length(ref))';
    fit_exp = fit(temp_x, double(ref), 'exp2');
    
    % Robust regression to relate reference to signal
    fitRef2 = robustfit(fit_exp(temp_x), theSig);
    
    % Create fitted signal
    fitSig = fit_exp(temp_x) * fitRef2(2) + fitRef2(1);
    
    % Correct signal
    F = theSig ./ fitSig;
    
    fprintf('Isosbestic correction complete\n');
else
    fprintf('Skipping isosbestic correction (sensor mode)\n');
    F = theSig;
end

% Remove bad periods
if ~isempty(rmStart)
    fprintf('Removing %d artifact periods...\n', length(rmStart));
    removeThis = false(size(theTime));
    for i = 1:length(rmStart)
        removeThis = removeThis | (theTime >= rmStart(i) & theTime <= rmEnd(i));
    end
    
    framesToRemove = sum(removeThis);
    if framesToRemove > 0 && framesToRemove < length(removeThis)
        F(removeThis) = [];
        theTime(removeThis) = [];
        fprintf('Removed %d frames (%.1f%%)\n', framesToRemove, ...
            100*framesToRemove/(framesToRemove + length(F)));
    else
        fprintf('WARNING: No frames removed (check rmStart/rmEnd values)\n');
    end
end

% Store for analysis
time = theTime;
clear theSig ref theTime;

fprintf('Final dataset: %d frames (%.1f minutes)\n', length(F), time(end)/60);

%% ========================================================================
%% BASELINE ESTIMATION (Slow-varying component)
%% ========================================================================

fprintf('\n========================================\n');
fprintf('ESTIMATING SLOW-VARYING BASELINE\n');
fprintf('========================================\n');

baseline_slow = estimateBaseline(F, time, BASELINE_WINDOW_MIN, BASELINE_PERCENTILE);

% Calculate dF/F relative to slow component
dF_F = (F - baseline_slow) ./ baseline_slow;

fprintf('Baseline estimation complete\n');
fprintf('Mean dF/F: %.4f\n', mean(dF_F));
fprintf('Std dF/F: %.4f\n', std(dF_F));

%% ========================================================================
%% Z-SCORE RELATIVE TO PRE-INJECTION BASELINE
%% ========================================================================

fprintf('\n========================================\n');
fprintf('Z-SCORING TO BASELINE PERIOD\n');
fprintf('========================================\n');

% Find pre-injection period
if injStart >= time(end)
    error('injStart (%.2f sec) is beyond end of data (%.2f sec)', injStart, time(end));
end

% Define baseline period while skipping early photobleaching
baselineStartSec = BASELINE_START_MIN * 60;
baselineMask = (time >= baselineStartSec) & (time < injStart);

if sum(baselineMask) < 100
    warning('Few baseline frames (%d). Consider lowering BASELINE_START_MIN or checking injStart.', sum(baselineMask));
end

baseline_dF_F = dF_F(baselineMask);

fprintf('Using baseline window: %.1f–%.1f min (%.0f–%.0f sec)\n', ...
    BASELINE_START_MIN, injStart/60, baselineStartSec, injStart);

% Calculate baseline statistics
baseline_stats.mean = mean(baseline_dF_F);
baseline_stats.std = std(baseline_dF_F);
baseline_stats.duration_min = sum(baselineMask) * mean(diff(time)) / 60;

% Z-score entire trace using baseline statistics
F_zscore = (dF_F - baseline_stats.mean) / baseline_stats.std;

fprintf('Baseline period: %.1f minutes (%d frames)\n', ...
    baseline_stats.duration_min, sum(baselineMask));
fprintf('Baseline mean dF/F: %.4f ± %.4f\n', ...
    baseline_stats.mean, baseline_stats.std);
fprintf('Post-baseline mean z-score: %.4f\n', mean(F_zscore(~baselineMask)));

%% ========================================================================
%% DETECT TRANSIENTS
%% ========================================================================

fprintf('\n========================================\n');
fprintf('DETECTING CALCIUM TRANSIENTS\n');
fprintf('========================================\n');

[transients, peak_times, peak_amps] = detectTransients(F_zscore, time, ...
    TRANSIENT_MIN_PROMINENCE, TRANSIENT_MIN_SEPARATION, TRANSIENT_MIN_WIDTH, ...
    TRANSIENT_RISE_MIN, TRANSIENT_RISE_MAX, TRANSIENT_DECAY_MIN, TRANSIENT_DECAY_MAX);

fprintf('Total transients detected: %d\n', length(transients.peak_times));
if ~isempty(transients.peak_amps)
    fprintf('Mean transient amplitude: %.2f z-score\n', mean(transients.peak_amps));
    fprintf('Transient frequency: %.2f events/min\n', ...
        length(transients.peak_times) / (time(end)/60));
else
    fprintf('No transients detected. Consider adjusting detection parameters.\n');
end

%% ========================================================================
%% QUANTIFY METRICS RELATIVE TO BASELINE
%% ========================================================================

fprintf('\n========================================\n');
fprintf('QUANTIFYING METRICS\n');
fprintf('========================================\n');

% Package time bins
timeBins.injStart = injStart;
timeBins.injEnd = injEnd;
timeBins.DSIStart = DSIStart;
timeBins.DSIEnd = DSIEnd;
timeBins.CageStart = CageStart;
timeBins.CageEnd = CageEnd;
timeBins.CagefoodStart = CagefoodStart;
timeBins.CagefoodEnd = CagefoodEnd;
timeBins.FoodStart = FoodStart;
timeBins.FoodEnd = FoodEnd;

% Comprehensive quantification
metrics = quantifyRelativeToBaseline(F_zscore, time, transients, timeBins,outFolder,BASELINE_START_MIN);

% Display summary
displayMetricsSummary(metrics);

%% ========================================================================
%% SAVE RESULTS - CSV FILES FOR GRAPHPAD PRISM
%% ========================================================================

fprintf('\n========================================\n');
fprintf('SAVING RESULTS\n');
fprintf('========================================\n');

% 1. Processed signals (downsampled and with period labels)
% Downsample by factor of 10
downsample_factor = PROCESSED_SIGNALS_DOWNSAMPLE;
downsample_idx = 1:downsample_factor:length(time);

time_ds = time(downsample_idx);
F_ds = F(downsample_idx);
baseline_slow_ds = baseline_slow(downsample_idx);
dF_F_ds = dF_F(downsample_idx);
F_zscore_ds = F_zscore(downsample_idx);

% Add period labels for each timepoint
period_labels_signals = cell(length(time_ds), 1);
for i = 1:length(time_ds)
    t = time_ds(i);
    
    % Determine which period this timepoint belongs to
    if t < timeBins.injStart
        period_labels_signals{i} = 'Baseline';
    elseif ~isempty(timeBins.injEnd) && t >= timeBins.injStart && t < timeBins.injEnd
        period_labels_signals{i} = 'Injection';
    elseif ~isempty(timeBins.injEnd) && ~isempty(timeBins.DSIStart) && t >= timeBins.injEnd && t < timeBins.DSIStart
        period_labels_signals{i} = 'Post_Injection';
    elseif ~isempty(timeBins.DSIStart) && ~isempty(timeBins.DSIEnd) && t >= timeBins.DSIStart && t < timeBins.DSIEnd
        period_labels_signals{i} = 'DSI';
    elseif ~isempty(timeBins.DSIEnd) && ~isempty(timeBins.CageStart) && t >= timeBins.DSIEnd && t < timeBins.CageStart
        period_labels_signals{i} = 'Post_DSI';
    elseif ~isempty(timeBins.CageStart) && ~isempty(timeBins.CageEnd) && t >= timeBins.CageStart && t < timeBins.CageEnd
        period_labels_signals{i} = 'Cage';
    elseif ~isempty(timeBins.CageEnd) && ~isempty(timeBins.CagefoodStart) && t >= timeBins.CageEnd && t < timeBins.CagefoodStart
        period_labels_signals{i} = 'ITI_1';
    elseif ~isempty(timeBins.CagefoodStart) && ~isempty(timeBins.CagefoodEnd) && t >= timeBins.CagefoodStart && t < timeBins.CagefoodEnd
        period_labels_signals{i} = 'Cage_Food';
    elseif ~isempty(timeBins.CagefoodEnd) && ~isempty(timeBins.FoodStart) && t >= timeBins.CagefoodEnd && t < timeBins.FoodStart
        period_labels_signals{i} = 'ITI_2';
    elseif ~isempty(timeBins.FoodStart) && t >= timeBins.FoodStart
        period_labels_signals{i} = 'Food';
    else
        period_labels_signals{i} = 'Unknown';
    end
end

% Create output table
output_signals = table(time_ds/60, period_labels_signals, F_ds, baseline_slow_ds, dF_F_ds, F_zscore_ds, ...
    'VariableNames', {'Time_min', 'Period', 'F_corrected', 'Baseline_slow', 'dF_F', 'F_zscore'});
writetable(output_signals, [outFolder 'processed_signals_' mouseName '.csv']);
fprintf('Saved: processed_signals_%s.csv (downsampled by %dx)\n', mouseName, downsample_factor);
% 2. Transient events
if ~isempty(transients.peak_times)
    output_transients = table(transients.peak_times/60, transients.peak_amps, ...
        transients.widths, transients.prominence, ...
        'VariableNames', {'Time_min', 'Amplitude_zscore', 'Width_sec', 'Prominence'});
    writetable(output_transients, [outFolder 'transients_' mouseName '.csv']);
    fprintf('Saved: transients_%s.csv\n', mouseName);
end


% 3. Binned time series (for plotting in Prism)
% Ensure all binned metrics are column vectors
fields = fieldnames(metrics.binned);
for i = 1:numel(fields)
    val = metrics.binned.(fields{i});
    if isrow(val)
        metrics.binned.(fields{i}) = val(:); % make column
    end
end

output_binned = table(metrics.binned.time_min, ...
    metrics.binned.period, ...
    metrics.binned.mean, ...
    metrics.binned.std, ...
    metrics.binned.sem, ...
    metrics.binned.delta_from_baseline, ...
    metrics.binned.transient_count, ...
    metrics.binned.auc_positive, ...
    metrics.binned.auc_negative, ...
    metrics.binned.auc_net, ...
    metrics.binned.auc_transient, ...
    metrics.binned.auc_sustained, ...
    'VariableNames', {'Time_Minutes', 'Period', 'Signal_Mean_Zscore', ...
    'Signal_Std_Zscore', 'Signal_SEM_Zscore', 'Delta_From_Baseline_Zscore', ...
    'Transient_Count_Per_Bin', ...
    'AUC_Positive_Per_Bin', 'AUC_Negative_Per_Bin', 'AUC_Net_Per_Bin', ...
    'AUC_Transient_Per_Bin', 'AUC_Sustained_Per_Bin'});
writetable(output_binned, [outFolder 'binned_1min_' mouseName '.csv']);
fprintf('Saved: binned_1min_%s.csv\n', mouseName);

% 4. Summary metrics by period
summaryTable = createSummaryMetricsTable(metrics);
writetable(summaryTable, [outFolder 'summary_metrics_' mouseName '.csv']);
fprintf('Saved: summary_metrics_%s.csv\n', mouseName);

% 5. Comparison metrics (deltas from baseline)
if isfield(metrics, 'comparison') && ~isempty(fieldnames(metrics.comparison))
    comparisonTable = createComparisonTable(metrics);
    writetable(comparisonTable, [outFolder 'comparison_to_baseline_' mouseName '.csv']);
    fprintf('Saved: comparison_to_baseline_%s.csv\n', mouseName);
end

% 6. Baseline period metrics
baselineTable = struct2table(metrics.baseline_period, 'AsArray', true);
writetable(baselineTable, [outFolder 'baseline_metrics_' mouseName '.csv']);
fprintf('Saved: baseline_metrics_%s.csv\n', mouseName);

% 7. Time-binned summaries (hourly and 30-min)
fprintf('\nGenerating time-binned summaries...\n');

% Hourly summary
createTimeBinnedSummary(F_zscore, time, transients, timeBins, baseline_stats, ...
    outFolder, mouseName, 60); % 60-minute bins

% 30-minute summary  
createTimeBinnedSummary(F_zscore, time, transients, timeBins, baseline_stats, ...
    outFolder, mouseName, 30); % 30-minute bins


% Save metrics structure (MATLAB format)
save([outFolder 'metrics_' mouseName '.mat'], 'metrics', 'baseline_stats', 'timeBins');

% Save summary report (text file)
saveSummaryReport(metrics, outFolder, mouseName);

fprintf('\nAll results saved to: %s\n', outFolder);

%% ========================================================================
%% VISUALIZATION
%% ========================================================================

%% ========================================================================
%% VISUALIZATION
%% ========================================================================

if PLOT_RESULTS
    fprintf('\n========================================\n');
    fprintf('GENERATING PLOTS\n');
    fprintf('========================================\n');
    
    visualizeResults(time, F, baseline_slow, F_zscore, transients, timeBins, ...
        metrics, outFolder, mouseName);
    
    % AUC component diagnostic
    fprintf('\nGenerating AUC diagnostic plots...\n');
    diagnostic_times = [20, 50, 105.5, 140]; % minutes - adjust as needed
    for i = 1:length(diagnostic_times)
        if diagnostic_times(i) < time(end)/60
            visualizeAUCComponents(F_zscore, time, transients, outFolder, mouseName, ...
                diagnostic_times(i), 3); % 3-minute window
        end
    end
    
    fprintf('Plots saved\n');
end

%% ========================================================================
%% FUNCTION DEFINITIONS
%% ========================================================================

function baseline = estimateBaseline(F, time, windowSize_min, percentile)
    % Estimate slow-varying baseline using percentile filtering
    
    fs = 1/mean(diff(time));
    windowSamples = round(windowSize_min * 60 * fs);
    
    baseline = zeros(size(F));
    halfWin = floor(windowSamples/2);
    
    for i = 1:length(F)
        startIdx = max(1, i - halfWin);
        endIdx = min(length(F), i + halfWin);
        
        % Use low percentile (resistant to upward transients)
        baseline(i) = prctile(F(startIdx:endIdx), percentile);
    end
    
    % Smooth the baseline estimate
    baseline = smoothdata(baseline, 'movmedian', round(windowSamples/10));
end

function [transients, peak_times, peak_amps] = detectTransients(F_zscore, time, ...
    minProminence, minPeakDistance_sec, minWidth_sec, ...
    riseMin, riseMax, decayMin, decayMax)
    % Detect calcium transients with GCaMP6s-appropriate kinetics
    
    fs = 1/mean(diff(time));
    minPeakDistance_samples = round(minPeakDistance_sec * fs);
    minWidth_samples = round(minWidth_sec * fs);
    
    % Find peaks
    [pks, locs, widths, proms] = findpeaks(F_zscore, ...
        'MinPeakProminence', minProminence, ...
        'MinPeakDistance', minPeakDistance_samples, ...
        'MinPeakWidth', minWidth_samples);
    
    % Validate kinetics
    validPeaks = false(size(locs));
    
    for i = 1:length(locs)
        peakIdx = locs(i);
        
        % Define search windows
        riseStart = max(1, peakIdx - round(2.0 * fs));
        decayEnd = min(length(F_zscore), peakIdx + round(3.0 * fs));
        
        % Find start of rise (crossing 0.5 SD)
        riseWindow = F_zscore(riseStart:peakIdx);
        riseStartIdx = find(riseWindow < 0.5, 1, 'last');
        if isempty(riseStartIdx)
            riseStartIdx = 1;
        end
        trueRiseStart = riseStart + riseStartIdx - 1;
        
        % Find end of decay (crossing 0.5 SD)
        decayWindow = F_zscore(peakIdx:decayEnd);
        decayEndIdx = find(decayWindow < 0.5, 1, 'first');
        if isempty(decayEndIdx)
            decayEndIdx = length(decayWindow);
        end
        trueDecayEnd = peakIdx + decayEndIdx - 1;
        
        % Calculate rise and decay times
        riseTime = time(peakIdx) - time(trueRiseStart);
        decayTime = time(trueDecayEnd) - time(peakIdx);
        
        % Validate kinetics for GCaMP6s
        validRise = riseTime >= riseMin && riseTime <= riseMax;
        validDecay = decayTime >= decayMin && decayTime <= decayMax;
        validRiseDecayRatio = decayTime > riseTime;
        
        validPeaks(i) = validRise && validDecay && validRiseDecayRatio;
    end
    
    % Filter to valid transients
    transients.locs = locs(validPeaks);
    transients.peak_times = time(locs(validPeaks));
    transients.peak_amps = pks(validPeaks);
    transients.widths = widths(validPeaks);
    transients.prominence = proms(validPeaks);
    
    peak_times = transients.peak_times;
    peak_amps = transients.peak_amps;
end

function metrics = quantifyRelativeToBaseline(F_zscore, time, transients, timeBins,outFolder,BASELINE_START_MIN)
    % Quantify all metrics relative to pre-injection baseline
    
    metrics.baseline_period = struct();
    metrics.comparison = struct();
    
    %% Baseline period (pre-injection)
    baselineStartSec = BASELINE_START_MIN * 60;
    injStart = timeBins.injStart;
    baselineMask = (time >= baselineStartSec) & (time < injStart);
    
    metrics.baseline_period.mean_zscore = mean(F_zscore(baselineMask));
    metrics.baseline_period.std_zscore = std(F_zscore(baselineMask));
    metrics.baseline_period.duration_min = sum(baselineMask) * mean(diff(time)) / 60;
    
    % Transients in baseline
    baseline_transient_mask = transients.peak_times < timeBins.injStart;
    baseline_transient_times = transients.peak_times(baseline_transient_mask);
    baseline_transient_amps = transients.peak_amps(baseline_transient_mask);
    
    metrics.baseline_period.transient_count = length(baseline_transient_times);
    metrics.baseline_period.transient_freq = length(baseline_transient_times) / metrics.baseline_period.duration_min;
    
    if ~isempty(baseline_transient_amps)
        metrics.baseline_period.transient_amp_mean = mean(baseline_transient_amps);
        metrics.baseline_period.transient_amp_median = median(baseline_transient_amps);
    else
        metrics.baseline_period.transient_amp_mean = NaN;
        metrics.baseline_period.transient_amp_median = NaN;
    end
    
    % ==================== BASELINE AUC CALCULATIONS ====================
    baseline_signal = F_zscore(baselineMask);
    baseline_time = time(baselineMask);
    
    % Global AUC (positive/negative split)
    positive_signal = max(baseline_signal, 0);
    negative_signal = abs(min(baseline_signal, 0));
    
    metrics.baseline_period.auc_positive_total = trapz(baseline_time, positive_signal);
    metrics.baseline_period.auc_negative_total = trapz(baseline_time, negative_signal);
    metrics.baseline_period.auc_net = trapz(baseline_time, baseline_signal);
    
    % Per-minute normalized
    metrics.baseline_period.auc_positive_per_min = metrics.baseline_period.auc_positive_total / metrics.baseline_period.duration_min;
    metrics.baseline_period.auc_negative_per_min = metrics.baseline_period.auc_negative_total / metrics.baseline_period.duration_min;
    metrics.baseline_period.auc_net_per_min = metrics.baseline_period.auc_net / metrics.baseline_period.duration_min;
    
    % Event-based transient AUC
    [transient_auc_total, transient_auc_mean] = calculateEventBasedAUC(...
        F_zscore, time, transients, baseline_transient_mask);
    
    metrics.baseline_period.auc_transient_total = transient_auc_total;
    metrics.baseline_period.auc_transient_per_event = transient_auc_mean;
    metrics.baseline_period.auc_transient_per_min = transient_auc_total / metrics.baseline_period.duration_min;
    
    % Sustained AUC (net minus transient)
    metrics.baseline_period.auc_sustained_net = metrics.baseline_period.auc_net - transient_auc_total;
    metrics.baseline_period.auc_sustained_per_min = metrics.baseline_period.auc_sustained_net / metrics.baseline_period.duration_min;


    %% Post-injection period
    if ~isempty(timeBins.injEnd) && timeBins.injEnd < max(time)
        if ~isempty(timeBins.DSIStart)
            postInjEnd = timeBins.DSIStart;
        elseif ~isempty(timeBins.FoodStart)
            postInjEnd = timeBins.FoodStart;
        else
            postInjEnd = max(time);
        end
        
        postInjMask = time >= timeBins.injEnd & time < postInjEnd;
        if sum(postInjMask) > 0
            metrics = compareToBaseline(metrics, F_zscore, time, transients, ...
                postInjMask, 'post_injection', metrics.baseline_period);
        end
    end
    
    %% Post-DSI period
    if ~isempty(timeBins.DSIEnd) && timeBins.DSIEnd < max(time)
        if ~isempty(timeBins.FoodStart)
            postDSIEnd = timeBins.FoodStart;
        else
            postDSIEnd = max(time);
        end
        
        postDSIMask = time >= timeBins.DSIEnd & time < postDSIEnd;
        if sum(postDSIMask) > 0
            metrics = compareToBaseline(metrics, F_zscore, time, transients, ...
                postDSIMask, 'post_DSI', metrics.baseline_period);
        end
    end
    
    %% Post-food period
    if ~isempty(timeBins.FoodStart) && timeBins.FoodStart < max(time)
        postFoodMask = time >= timeBins.FoodStart;
        if sum(postFoodMask) > 0
            metrics = compareToBaseline(metrics, F_zscore, time, transients, ...
                postFoodMask, 'post_food', metrics.baseline_period);
        end
    end
    
    %% Binned time series with period identification and transient metrics
[binned, bin_centers] = binSignal(F_zscore, time, 1.0);

% Create period identifier for each bin
period_labels = cell(length(bin_centers), 1);
transient_counts_per_bin = zeros(length(bin_centers), 1);

% Initialize AUC arrays
auc_positive_per_bin = zeros(length(bin_centers), 1);
auc_negative_per_bin = zeros(length(bin_centers), 1);
auc_net_per_bin = zeros(length(bin_centers), 1);
auc_transient_per_bin = zeros(length(bin_centers), 1);
auc_sustained_per_bin = zeros(length(bin_centers), 1);

for i = 1:length(bin_centers)
    bin_time = bin_centers(i);
    
    % Determine which period this bin belongs to
    if bin_time < timeBins.injStart
        period_labels{i} = 'Baseline';
    elseif ~isempty(timeBins.injEnd) && bin_time >= timeBins.injStart && bin_time < timeBins.injEnd
        period_labels{i} = 'Injection';
    elseif ~isempty(timeBins.injEnd) && ~isempty(timeBins.DSIStart) && bin_time >= timeBins.injEnd && bin_time < timeBins.DSIStart
        period_labels{i} = 'Post_Injection';
    elseif ~isempty(timeBins.DSIStart) && ~isempty(timeBins.DSIEnd) && bin_time >= timeBins.DSIStart && bin_time < timeBins.DSIEnd
        period_labels{i} = 'DSI';
    elseif ~isempty(timeBins.DSIEnd) && ~isempty(timeBins.CageStart) && bin_time >= timeBins.DSIEnd && bin_time < timeBins.CageStart
        period_labels{i} = 'Post_DSI';
    elseif ~isempty(timeBins.CageStart) && ~isempty(timeBins.CageEnd) && bin_time >= timeBins.CageStart && bin_time < timeBins.CageEnd
        period_labels{i} = 'Cage';
    elseif ~isempty(timeBins.CageEnd) && ~isempty(timeBins.CagefoodStart) && bin_time >= timeBins.CageEnd && bin_time < timeBins.CagefoodStart
        period_labels{i} = 'ITI_1';
    elseif ~isempty(timeBins.CagefoodStart) && ~isempty(timeBins.CagefoodEnd) && bin_time >= timeBins.CagefoodStart && bin_time < timeBins.CagefoodEnd
        period_labels{i} = 'Cage_Food';
    elseif ~isempty(timeBins.CagefoodEnd) && ~isempty(timeBins.FoodStart) && bin_time >= timeBins.CagefoodEnd && bin_time < timeBins.FoodStart
        period_labels{i} = 'ITI_2';
    elseif ~isempty(timeBins.FoodStart) && bin_time >= timeBins.FoodStart
        period_labels{i} = 'Food';
    else
        period_labels{i} = 'Unknown';
    end
    
    % Count transients in this bin
    bin_start = bin_centers(i) - 30;  % 30 sec before center (for 1-min bins)
    bin_end = bin_centers(i) + 30;    % 30 sec after center
    transients_in_bin_mask = transients.peak_times >= bin_start & transients.peak_times < bin_end;
    transient_counts_per_bin(i) = sum(transients_in_bin_mask);
    
    % Calculate ALL AUC metrics for this bin
    time_mask = time >= bin_start & time < bin_end;
    if sum(time_mask) > 0
        bin_signal = F_zscore(time_mask);
        bin_time_vals = time(time_mask);
        
        % Global AUC (positive/negative split)
        positive_signal = max(bin_signal, 0);
        negative_signal = abs(min(bin_signal, 0));
        
        auc_positive_per_bin(i) = trapz(bin_time_vals, positive_signal);
        auc_negative_per_bin(i) = trapz(bin_time_vals, negative_signal);
        auc_net_per_bin(i) = trapz(bin_time_vals, bin_signal);
        
        % Event-based transient AUC for this bin
        [transient_auc_bin, ~] = calculateEventBasedAUC(...
            F_zscore, time, transients, transients_in_bin_mask);
        auc_transient_per_bin(i) = transient_auc_bin;
        
        % Sustained AUC (net minus transient)
        auc_sustained_per_bin(i) = auc_net_per_bin(i) - transient_auc_bin;
    else
        auc_positive_per_bin(i) = 0;
        auc_negative_per_bin(i) = 0;
        auc_net_per_bin(i) = 0;
        auc_transient_per_bin(i) = 0;
        auc_sustained_per_bin(i) = 0;
    end
end

% ----------------- Diagnostic & save figure (around a chosen bin) -----------------
% Adjust center_min as needed; 
diagnostic_center_min = 105.5;
diagnostic_center_sec = diagnostic_center_min * 60;
win = 90;  % seconds around center to show

time_mask_diag = time >= (diagnostic_center_sec - win) & time <= (diagnostic_center_sec + win);
if sum(time_mask_diag) > 10
    fig = figure('Visible','off');
    plot(time(time_mask_diag)/60, F_zscore(time_mask_diag), 'k'); hold on;
    xlabel('Time (min)'); ylabel('Z-score');
    title(sprintf('Diagnostic around %.1f min', diagnostic_center_min));

    % plot global-detected peaks (if any)
    if isfield(transients,'peak_times') && ~isempty(transients.peak_times)
        pk_mask = transients.peak_times >= (diagnostic_center_sec - win) & transients.peak_times <= (diagnostic_center_sec + win);
        if any(pk_mask)
            plot(transients.peak_times(pk_mask)/60, transients.peak_amps(pk_mask), 'ro','MarkerFaceColor','r');
        end
    end

    % bin boundaries and thresholds for reference
    xline(diagnostic_center_min - 0.5, 'b--');
    xline(diagnostic_center_min + 0.5, 'b--');
    yline(2, 'r--', '2 SD');
    yline(1, 'm--', '1 SD');

    % build legend robustly
    legend_entries = {'Z-scored signal'};
    if exist('pk_mask','var') && any(pk_mask)
        legend_entries{end+1} = 'Detected peaks';
    end
    legend(legend_entries, 'Location', 'best');

    % Save into the output folder provided by the main script
    try
        saveas(fig, fullfile(outFolder, sprintf('diagnostic_transients_around_%.1fmin.png', diagnostic_center_min)));
    catch ME
        % Fallback: save to current folder if saving to outFolder fails
        warning('FiberPhotometry:SaveError', 'Could not save to outFolder: %s. Saving to pwd instead.', ME.message);
        saveas(fig, sprintf('diagnostic_transients_around_%.1fmin.png', diagnostic_center_min));
    end
    close(fig);
end

metrics.binned.time_min = bin_centers / 60;
metrics.binned.mean = binned.mean;
metrics.binned.std = binned.std;
metrics.binned.sem = binned.sem;
metrics.binned.delta_from_baseline = binned.mean - metrics.baseline_period.mean_zscore;
metrics.binned.period = period_labels;
metrics.binned.transient_count = transient_counts_per_bin;

% AUC metrics per bin
metrics.binned.auc_positive = auc_positive_per_bin;
metrics.binned.auc_negative = auc_negative_per_bin;
metrics.binned.auc_net = auc_net_per_bin;
metrics.binned.auc_transient = auc_transient_per_bin;
metrics.binned.auc_sustained = auc_sustained_per_bin;

end

function metrics = compareToBaseline(metrics, F_zscore, time, transients, periodMask, periodName, baseline_ref)
    % Compare a given period to the baseline reference
    
    period_data = F_zscore(periodMask);
    period_time = time(periodMask);
    
    if isempty(period_data)
        return;
    end
    
    duration_min = length(period_data) * mean(diff(time)) / 60;
    
    %% Mean signal level
    metrics.comparison.(periodName).mean_zscore = mean(period_data);
    metrics.comparison.(periodName).std_zscore = std(period_data);
    metrics.comparison.(periodName).duration_min = duration_min;
    
    % Delta from baseline
    metrics.comparison.(periodName).delta_mean = ...
        metrics.comparison.(periodName).mean_zscore - baseline_ref.mean_zscore;
    
    % Effect size (Cohen's d)
    pooled_std = sqrt((baseline_ref.std_zscore^2 + metrics.comparison.(periodName).std_zscore^2) / 2);
    metrics.comparison.(periodName).cohens_d = ...
        metrics.comparison.(periodName).delta_mean / pooled_std;
    
    %% Transient metrics
    transient_mask = transients.peak_times >= min(period_time) & ...
                     transients.peak_times <= max(period_time);
    period_transient_times = transients.peak_times(transient_mask);
    period_transient_amps = transients.peak_amps(transient_mask);
    
    metrics.comparison.(periodName).transient_count = length(period_transient_times);
    metrics.comparison.(periodName).transient_freq = length(period_transient_times) / duration_min;
    
    % Delta from baseline
    metrics.comparison.(periodName).delta_transient_freq = ...
        metrics.comparison.(periodName).transient_freq - baseline_ref.transient_freq;
    
    if baseline_ref.transient_freq > 0
        metrics.comparison.(periodName).delta_transient_freq_pct = ...
            (metrics.comparison.(periodName).delta_transient_freq / baseline_ref.transient_freq) * 100;
    else
        metrics.comparison.(periodName).delta_transient_freq_pct = NaN;
    end
    
    % Amplitude changes
    if ~isempty(period_transient_amps)
        metrics.comparison.(periodName).transient_amp_mean = mean(period_transient_amps);
        metrics.comparison.(periodName).delta_transient_amp = ...
            mean(period_transient_amps) - baseline_ref.transient_amp_mean;
    else
        metrics.comparison.(periodName).transient_amp_mean = NaN;
        metrics.comparison.(periodName).delta_transient_amp = NaN;
    end
    
    %% AUC metrics
    % Global AUC (positive/negative split)
    positive_signal = max(period_data, 0);
    negative_signal = abs(min(period_data, 0));
    
    metrics.comparison.(periodName).auc_positive_total = trapz(period_time, positive_signal);
    metrics.comparison.(periodName).auc_negative_total = trapz(period_time, negative_signal);
    metrics.comparison.(periodName).auc_net = trapz(period_time, period_data);
    
    % Per-minute normalized
    metrics.comparison.(periodName).auc_positive_per_min = metrics.comparison.(periodName).auc_positive_total / duration_min;
    metrics.comparison.(periodName).auc_negative_per_min = metrics.comparison.(periodName).auc_negative_total / duration_min;
    metrics.comparison.(periodName).auc_net_per_min = metrics.comparison.(periodName).auc_net / duration_min;
    
    % Event-based transient AUC
    [transient_auc_total, transient_auc_mean] = calculateEventBasedAUC(...
        F_zscore, time, transients, transient_mask);
    
    metrics.comparison.(periodName).auc_transient_total = transient_auc_total;
    metrics.comparison.(periodName).auc_transient_per_event = transient_auc_mean;
    metrics.comparison.(periodName).auc_transient_per_min = transient_auc_total / duration_min;
    
    % Sustained AUC
    metrics.comparison.(periodName).auc_sustained_net = metrics.comparison.(periodName).auc_net - transient_auc_total;
    metrics.comparison.(periodName).auc_sustained_per_min = metrics.comparison.(periodName).auc_sustained_net / duration_min;
    
    % ===== DELTAS FROM BASELINE =====
    % Positive AUC delta
    metrics.comparison.(periodName).delta_auc_positive_per_min = ...
        metrics.comparison.(periodName).auc_positive_per_min - baseline_ref.auc_positive_per_min;
    
    % Negative AUC delta
    metrics.comparison.(periodName).delta_auc_negative_per_min = ...
        metrics.comparison.(periodName).auc_negative_per_min - baseline_ref.auc_negative_per_min;
    
    % Net AUC delta
    metrics.comparison.(periodName).delta_auc_net_per_min = ...
        metrics.comparison.(periodName).auc_net_per_min - baseline_ref.auc_net_per_min;
    
    % Transient AUC delta
    metrics.comparison.(periodName).delta_auc_transient_per_min = ...
        metrics.comparison.(periodName).auc_transient_per_min - baseline_ref.auc_transient_per_min;
    
    % Sustained AUC delta
    metrics.comparison.(periodName).delta_auc_sustained_per_min = ...
        metrics.comparison.(periodName).auc_sustained_per_min - baseline_ref.auc_sustained_per_min;
    
    % Percent changes (where applicable)
    if baseline_ref.auc_positive_per_min > 0
        metrics.comparison.(periodName).delta_auc_positive_pct = ...
            (metrics.comparison.(periodName).delta_auc_positive_per_min / baseline_ref.auc_positive_per_min) * 100;
    else
        metrics.comparison.(periodName).delta_auc_positive_pct = NaN;
    end
    
    if baseline_ref.auc_net_per_min ~= 0
        metrics.comparison.(periodName).delta_auc_net_pct = ...
            (metrics.comparison.(periodName).delta_auc_net_per_min / baseline_ref.auc_net_per_min) * 100;
    else
        metrics.comparison.(periodName).delta_auc_net_pct = NaN;
    end
    
    %% Time above threshold
    metrics.comparison.(periodName).pct_time_above_2SD = ...
        sum(period_data > 2.0) / length(period_data) * 100;
    
    % Expected ~2.3% for normal distribution
    metrics.comparison.(periodName).delta_pct_above_2SD = ...
        metrics.comparison.(periodName).pct_time_above_2SD - 2.3;
end

function [binned, bin_centers] = binSignal(signal, time, binSize_min)
    % Bin signal into equal time bins
    
    binSize_sec = binSize_min * 60;
    
    bin_edges = min(time):binSize_sec:max(time);
    if bin_edges(end) < max(time)
        bin_edges = [bin_edges, max(time)];
    end
    
    nBins = length(bin_edges) - 1;
    
    binned.mean = zeros(1, nBins);
    binned.std = zeros(1, nBins);
    binned.sem = zeros(1, nBins);
    binned.n_samples = zeros(1, nBins);
    bin_centers = zeros(1, nBins);
    
    for i = 1:nBins
        mask = time >= bin_edges(i) & time < bin_edges(i+1);
        bin_data = signal(mask);
        
        if ~isempty(bin_data)
            binned.mean(i) = mean(bin_data);
            binned.std(i) = std(bin_data);
            binned.sem(i) = std(bin_data) / sqrt(length(bin_data));
            binned.n_samples(i) = length(bin_data);
        end
        
        bin_centers(i) = (bin_edges(i) + bin_edges(i+1)) / 2;
    end
end

function summaryTable = createSummaryMetricsTable(metrics)
    % Create table with summary metrics for each period
    
    % Start with baseline
    periods = {'Baseline'};
    mean_zscore = metrics.baseline_period.mean_zscore;
    transient_freq = metrics.baseline_period.transient_freq;
    transient_amp = metrics.baseline_period.transient_amp_mean;
    
    % AUC metrics
    auc_positive_per_min = metrics.baseline_period.auc_positive_per_min;
    auc_negative_per_min = metrics.baseline_period.auc_negative_per_min;
    auc_net_per_min = metrics.baseline_period.auc_net_per_min;
    auc_transient_per_min = metrics.baseline_period.auc_transient_per_min;
    auc_sustained_per_min = metrics.baseline_period.auc_sustained_per_min;
    
    % Add comparison periods
    if isfield(metrics, 'comparison')
        compPeriods = fieldnames(metrics.comparison);
        for i = 1:length(compPeriods)
            periods{end+1} = compPeriods{i};
            mean_zscore(end+1) = metrics.comparison.(compPeriods{i}).mean_zscore;
            transient_freq(end+1) = metrics.comparison.(compPeriods{i}).transient_freq;
            transient_amp(end+1) = metrics.comparison.(compPeriods{i}).transient_amp_mean;
            
            % AUC metrics
            auc_positive_per_min(end+1) = metrics.comparison.(compPeriods{i}).auc_positive_per_min;
            auc_negative_per_min(end+1) = metrics.comparison.(compPeriods{i}).auc_negative_per_min;
            auc_net_per_min(end+1) = metrics.comparison.(compPeriods{i}).auc_net_per_min;
            auc_transient_per_min(end+1) = metrics.comparison.(compPeriods{i}).auc_transient_per_min;
            auc_sustained_per_min(end+1) = metrics.comparison.(compPeriods{i}).auc_sustained_per_min;
        end
    end
    
    summaryTable = table(periods', mean_zscore', transient_freq', transient_amp', ...
        auc_positive_per_min', auc_negative_per_min', auc_net_per_min', ...
        auc_transient_per_min', auc_sustained_per_min', ...
        'VariableNames', {'Period', 'Mean_zscore', 'Transient_freq', 'Transient_amplitude', ...
        'AUC_positive_per_min', 'AUC_negative_per_min', 'AUC_net_per_min', ...
        'AUC_transient_per_min', 'AUC_sustained_per_min'});
end

function comparisonTable = createComparisonTable(metrics)
    % Create table with delta metrics (changes from baseline)
    
    periods = {};
    delta_mean = [];
    cohens_d = [];
    delta_freq = [];
    delta_freq_pct = [];
    delta_amp = [];
    
    % AUC deltas
    delta_auc_positive = [];
    delta_auc_negative = [];
    delta_auc_net = [];
    delta_auc_net_pct = [];
    delta_auc_transient = [];
    delta_auc_sustained = [];
    
    compPeriods = fieldnames(metrics.comparison);
    for i = 1:length(compPeriods)
        periods{end+1} = compPeriods{i};
        delta_mean(end+1) = metrics.comparison.(compPeriods{i}).delta_mean;
        cohens_d(end+1) = metrics.comparison.(compPeriods{i}).cohens_d;
        delta_freq(end+1) = metrics.comparison.(compPeriods{i}).delta_transient_freq;
        delta_freq_pct(end+1) = metrics.comparison.(compPeriods{i}).delta_transient_freq_pct;
        delta_amp(end+1) = metrics.comparison.(compPeriods{i}).delta_transient_amp;
        
        % AUC deltas
        delta_auc_positive(end+1) = metrics.comparison.(compPeriods{i}).delta_auc_positive_per_min;
        delta_auc_negative(end+1) = metrics.comparison.(compPeriods{i}).delta_auc_negative_per_min;
        delta_auc_net(end+1) = metrics.comparison.(compPeriods{i}).delta_auc_net_per_min;
        delta_auc_net_pct(end+1) = metrics.comparison.(compPeriods{i}).delta_auc_net_pct;
        delta_auc_transient(end+1) = metrics.comparison.(compPeriods{i}).delta_auc_transient_per_min;
        delta_auc_sustained(end+1) = metrics.comparison.(compPeriods{i}).delta_auc_sustained_per_min;
    end
    
    comparisonTable = table(periods', delta_mean', cohens_d', delta_freq', ...
        delta_freq_pct', delta_amp', ...
        delta_auc_positive', delta_auc_negative', delta_auc_net', delta_auc_net_pct', ...
        delta_auc_transient', delta_auc_sustained', ...
        'VariableNames', {'Period', 'Delta_mean_zscore', 'Cohens_d', ...
        'Delta_transient_freq', 'Delta_freq_percent', 'Delta_transient_amp', ...
        'Delta_AUC_positive_per_min', 'Delta_AUC_negative_per_min', ...
        'Delta_AUC_net_per_min', 'Delta_AUC_net_percent', ...
        'Delta_AUC_transient_per_min', 'Delta_AUC_sustained_per_min'});
end

function displayMetricsSummary(metrics)
    % Display summary of key metrics
    
    fprintf('\n--- BASELINE PERIOD ---\n');
    fprintf('Duration: %.1f min\n', metrics.baseline_period.duration_min);
    fprintf('Mean z-score: %.3f\n', metrics.baseline_period.mean_zscore);
    fprintf('Transient count: %d\n', metrics.baseline_period.transient_count);
    fprintf('Transient frequency: %.2f events/min\n', metrics.baseline_period.transient_freq);
    fprintf('Mean transient amplitude: %.2f z-score\n', metrics.baseline_period.transient_amp_mean);
    
    % AUC metrics
    fprintf('AUC positive per min: %.2f\n', metrics.baseline_period.auc_positive_per_min);
    fprintf('AUC negative per min: %.2f\n', metrics.baseline_period.auc_negative_per_min);
    fprintf('AUC net per min: %.2f\n', metrics.baseline_period.auc_net_per_min);
    fprintf('AUC transient per min: %.2f\n', metrics.baseline_period.auc_transient_per_min);
    fprintf('AUC sustained per min: %.2f\n', metrics.baseline_period.auc_sustained_per_min);
    
    % Display each comparison period
    if isfield(metrics, 'comparison')
        periods = fieldnames(metrics.comparison);
        
        for i = 1:length(periods)
            period = periods{i};
            fprintf('\n--- %s ---\n', upper(strrep(period, '_', ' ')));
            fprintf('Duration: %.1f min\n', metrics.comparison.(period).duration_min);
            fprintf('Mean z-score: %.3f (Δ = %+.3f, d = %.2f)\n', ...
                metrics.comparison.(period).mean_zscore, ...
                metrics.comparison.(period).delta_mean, ...
                metrics.comparison.(period).cohens_d);
            fprintf('Transient frequency: %.2f events/min (Δ = %+.2f, %+.1f%%)\n', ...
                metrics.comparison.(period).transient_freq, ...
                metrics.comparison.(period).delta_transient_freq, ...
                metrics.comparison.(period).delta_transient_freq_pct);
            fprintf('Mean transient amplitude: %.2f z-score (Δ = %+.2f)\n', ...
                metrics.comparison.(period).transient_amp_mean, ...
                metrics.comparison.(period).delta_transient_amp);
            
            % AUC metrics
            fprintf('AUC positive per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_positive_per_min, ...
                metrics.comparison.(period).delta_auc_positive_per_min);
            fprintf('AUC negative per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_negative_per_min, ...
                metrics.comparison.(period).delta_auc_negative_per_min);
            fprintf('AUC net per min: %.2f (Δ = %+.2f, %+.1f%%)\n', ...
                metrics.comparison.(period).auc_net_per_min, ...
                metrics.comparison.(period).delta_auc_net_per_min, ...
                metrics.comparison.(period).delta_auc_net_pct);
            fprintf('AUC transient per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_transient_per_min, ...
                metrics.comparison.(period).delta_auc_transient_per_min);
            fprintf('AUC sustained per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_sustained_per_min, ...
                metrics.comparison.(period).delta_auc_sustained_per_min);
            
            fprintf('Time above 2SD: %.1f%% (Δ = %+.1f%%)\n', ...
                metrics.comparison.(period).pct_time_above_2SD, ...
                metrics.comparison.(period).delta_pct_above_2SD);
        end
    end
end

function saveSummaryReport(metrics, outFolder, mouseName)
    % Save text summary report
    
    fid = fopen([outFolder 'summary_report_' mouseName '.txt'], 'w');
    
    fprintf(fid, '========================================\n');
    fprintf(fid, 'PHOTOMETRY ANALYSIS SUMMARY\n');
    fprintf(fid, 'Mouse: %s\n', mouseName);
    fprintf(fid, 'Date: %s\n', datestr(now));
    fprintf(fid, '========================================\n\n');
    
    fprintf(fid, '--- BASELINE PERIOD ---\n');
    fprintf(fid, 'Duration: %.1f min\n', metrics.baseline_period.duration_min);
    fprintf(fid, 'Mean z-score: %.3f ± %.3f\n', metrics.baseline_period.mean_zscore, metrics.baseline_period.std_zscore);
    fprintf(fid, 'Transient count: %d\n', metrics.baseline_period.transient_count);
    fprintf(fid, 'Transient frequency: %.2f events/min\n', metrics.baseline_period.transient_freq);
    fprintf(fid, 'Mean transient amplitude: %.2f z-score\n', metrics.baseline_period.transient_amp_mean);
    
    % AUC metrics
    fprintf(fid, 'AUC positive per min: %.2f\n', metrics.baseline_period.auc_positive_per_min);
    fprintf(fid, 'AUC negative per min: %.2f\n', metrics.baseline_period.auc_negative_per_min);
    fprintf(fid, 'AUC net per min: %.2f\n', metrics.baseline_period.auc_net_per_min);
    fprintf(fid, 'AUC transient per min: %.2f\n', metrics.baseline_period.auc_transient_per_min);
    fprintf(fid, 'AUC sustained per min: %.2f\n\n', metrics.baseline_period.auc_sustained_per_min);
    
    if isfield(metrics, 'comparison')
        periods = fieldnames(metrics.comparison);
        for i = 1:length(periods)
            period = periods{i};
            fprintf(fid, '--- %s ---\n', upper(strrep(period, '_', ' ')));
            fprintf(fid, 'Duration: %.1f min\n', metrics.comparison.(period).duration_min);
            fprintf(fid, 'Mean z-score: %.3f (Δ = %+.3f, Cohen''s d = %.2f)\n', ...
                metrics.comparison.(period).mean_zscore, ...
                metrics.comparison.(period).delta_mean, ...
                metrics.comparison.(period).cohens_d);
            fprintf(fid, 'Transient frequency: %.2f events/min (Δ = %+.2f, %+.1f%%)\n', ...
                metrics.comparison.(period).transient_freq, ...
                metrics.comparison.(period).delta_transient_freq, ...
                metrics.comparison.(period).delta_transient_freq_pct);
            fprintf(fid, 'Mean transient amplitude: %.2f z-score (Δ = %+.2f)\n', ...
                metrics.comparison.(period).transient_amp_mean, ...
                metrics.comparison.(period).delta_transient_amp);
            
            % AUC metrics
            fprintf(fid, 'AUC positive per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_positive_per_min, ...
                metrics.comparison.(period).delta_auc_positive_per_min);
            fprintf(fid, 'AUC negative per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_negative_per_min, ...
                metrics.comparison.(period).delta_auc_negative_per_min);
            fprintf(fid, 'AUC net per min: %.2f (Δ = %+.2f, %+.1f%%)\n', ...
                metrics.comparison.(period).auc_net_per_min, ...
                metrics.comparison.(period).delta_auc_net_per_min, ...
                metrics.comparison.(period).delta_auc_net_pct);
            fprintf(fid, 'AUC transient per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_transient_per_min, ...
                metrics.comparison.(period).delta_auc_transient_per_min);
            fprintf(fid, 'AUC sustained per min: %.2f (Δ = %+.2f)\n', ...
                metrics.comparison.(period).auc_sustained_per_min, ...
                metrics.comparison.(period).delta_auc_sustained_per_min);
            
            fprintf(fid, 'Time above 2SD: %.1f%% (Δ = %+.1f%%)\n\n', ...
                metrics.comparison.(period).pct_time_above_2SD, ...
                metrics.comparison.(period).delta_pct_above_2SD);
        end
    end
    
    fclose(fid);
end

function visualizeResults(time, F, baseline_slow, F_zscore, transients, timeBins, ...
    metrics, outFolder, mouseName)
    % Create comprehensive visualization
    
    time_min = time / 60;
    
    % Figure 1: Signal processing steps
    figure('Position', [100 100 1200 800]);
    
    subplot(4,1,1);
    plot(time_min, F, 'k', 'LineWidth', 0.5);
    hold on;
    plot(time_min, baseline_slow, 'r', 'LineWidth', 1.5);
    ylabel('Fluorescence (a.u.)');
    title(['Signal Processing - ' mouseName]);
    legend('Corrected F', 'Slow baseline', 'Location', 'best');
    addEventLines(timeBins, gca);
    grid on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:1:ceil(max(time_min));
    
    subplot(4,1,2);
    dF_F = (F - baseline_slow) ./ baseline_slow;
    plot(time_min, dF_F * 100, 'k', 'LineWidth', 0.5);
    ylabel('\DeltaF/F (%)');
    legend('dF/F', 'Location', 'best');
    addEventLines(timeBins, gca);
    grid on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:1:ceil(max(time_min));
    
    subplot(4,1,3);
    plot(time_min, F_zscore, 'k', 'LineWidth', 0.5);
    hold on;
    yline(0, 'b--', 'Baseline mean');
    yline(2, 'r--', '2 SD');
    yline(-2, 'r--', '-2 SD');
    ylabel('Z-score');
    legend('Z-scored signal', 'Location', 'best');
    addEventLines(timeBins, gca);
    grid on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:1:ceil(max(time_min));
    
    subplot(4,1,4);
    plot(time_min, F_zscore, 'k', 'LineWidth', 0.5);
    hold on;
    if ~isempty(transients.peak_times)
        plot(transients.peak_times/60, transients.peak_amps, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    end
    yline(0, 'b--', 'Baseline');
    ylabel('Z-score');
    xlabel('Time (min)');
    legend('Signal', 'Detected transients', 'Location', 'best');
    addEventLines(timeBins, gca);
    grid on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:1:ceil(max(time_min));
    
    saveas(gcf, [outFolder 'signal_processing_' mouseName '.png']);
    
    % Figure 2: Binned data with statistics
    figure('Position', [100 100 1200 600]);
    
    subplot(2,1,1);
    errorbar(metrics.binned.time_min, metrics.binned.mean, metrics.binned.sem, 'k.-');
    hold on;
    yline(0, 'b--', 'Baseline', 'LineWidth', 1.5);
    ylabel('Mean Z-score');
    title(['1-Minute Binned Signal - ' mouseName]);
    addEventLines(timeBins, gca);
    grid on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:1:ceil(max(metrics.binned.time_min));
    
    subplot(2,1,2);
    plot(metrics.binned.time_min, metrics.binned.delta_from_baseline, 'k.-');
    hold on;
    yline(0, 'r--', 'Zero change', 'LineWidth', 1.5);
    ylabel('\Delta from baseline (z-score)');
    xlabel('Time (min)');
    addEventLines(timeBins, gca);
    grid on;
    ax = gca;
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:1:ceil(max(metrics.binned.time_min));
    
    saveas(gcf, [outFolder 'binned_analysis_' mouseName '.png']);
    
    close all;
end

function addEventLines(timeBins, ax)
    % Add vertical lines for experimental events
    
    ylims = ylim(ax);
    
    if ~isempty(timeBins.injStart)
        line([timeBins.injStart timeBins.injStart]/60, ylims, ...
            'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1);
        text(timeBins.injStart/60, ylims(2)*0.95, 'Inj', 'FontSize', 8);
    end
    
    if ~isempty(timeBins.DSIStart)
        line([timeBins.DSIStart timeBins.DSIStart]/60, ylims, ...
            'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1);
        text(timeBins.DSIStart/60, ylims(2)*0.95, 'DSI', 'FontSize', 8);
    end
    
    if ~isempty(timeBins.FoodStart)
        line([timeBins.FoodStart timeBins.FoodStart]/60, ylims, ...
            'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1);
        text(timeBins.FoodStart/60, ylims(2)*0.95, 'Food', 'FontSize', 8);
    end
end

function visualizeAUCComponents(F_zscore, time, transients, outFolder, mouseName, diagnostic_center_min, window_min)
    % Visualize different AUC components for verification
    % 
    % Creates a figure showing:
    % - Global positive/negative AUC (shaded areas above/below baseline)
    % - Event-based transient AUC (with local baselines)
    % - Individual transient boundaries
    
    if nargin < 7
        window_min = 3; % default 3-minute window
    end
    
    % Convert to time indices
    diagnostic_center_sec = diagnostic_center_min * 60;
    window_sec = window_min * 60 / 2; % half window on each side
    
    % Find time window
    time_mask = time >= (diagnostic_center_sec - window_sec) & ...
                time <= (diagnostic_center_sec + window_sec);
    
    if sum(time_mask) < 10
        warning('Insufficient data points in diagnostic window');
        return;
    end
    
    time_window = time(time_mask);
    signal_window = F_zscore(time_mask);
    
    % Find transients in this window
    transient_mask = transients.peak_times >= (diagnostic_center_sec - window_sec) & ...
                     transients.peak_times <= (diagnostic_center_sec + window_sec);
    transients_in_window = structfun(@(x) x(transient_mask), transients, 'UniformOutput', false);
    
    % Create figure
    fig = figure('Position', [100, 100, 1400, 900], 'Visible', 'off');
    
    %% SUBPLOT 1: Global AUC (Positive/Negative split)
    subplot(3,1,1);
    hold on;
    
    % Plot baseline
    plot(time_window/60, zeros(size(time_window)), 'k--', 'LineWidth', 1.5);
    
    % Shade positive AUC (above baseline)
    pos_signal = max(signal_window, 0);
    area(time_window/60, pos_signal, 'FaceColor', [0.2 0.8 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Shade negative AUC (below baseline)
    neg_signal = min(signal_window, 0);
    area(time_window/60, neg_signal, 'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Plot signal on top
    plot(time_window/60, signal_window, 'k-', 'LineWidth', 1.5);
    
    % Calculate AUC values for display
    auc_pos = trapz(time_window, pos_signal);
    auc_neg = trapz(time_window, abs(neg_signal));
    auc_net = trapz(time_window, signal_window);
    
    ylabel('Z-score');
    title(sprintf('Global AUC: Positive (green) vs Negative (red)\nAUC_{pos}=%.2f, AUC_{neg}=%.2f, AUC_{net}=%.2f', ...
        auc_pos, auc_neg, auc_net), 'FontSize', 11);
    legend({'Baseline (z=0)', 'Positive AUC', 'Negative AUC', 'Signal'}, 'Location', 'best');
    grid on;
    xlim([time_window(1)/60, time_window(end)/60]);
    ax = gca;
    ax.XMinorTick = 'on';
    x_start = floor(time_window(1)/60);
    x_end = ceil(time_window(end)/60);
    ax.XAxis.MinorTickValues = x_start:1:x_end;
    
    %% SUBPLOT 2: Event-Based Transient AUC
    subplot(3,1,2);
    hold on;
    
    % Plot baseline
    plot(time_window/60, zeros(size(time_window)), 'k--', 'LineWidth', 1.5);
    
    % Plot signal
    plot(time_window/60, signal_window, 'k-', 'LineWidth', 1.5);
    
    % Process each transient
    total_transient_auc = 0;
    colors = lines(length(transients_in_window.peak_times));
    
    for i = 1:length(transients_in_window.peak_times)
        peak_time = transients_in_window.peak_times(i);
        peak_idx_global = transients_in_window.locs(i);
        
        % Define transient window (±3 sec or based on width)
        transient_width = transients_in_window.widths(i);
        window_size = max(3.0, transient_width * 2); % at least 6 sec total
        
        t_start = peak_time - window_size;
        t_end = peak_time + window_size;
        
        % Find indices in global time array
        event_mask_global = time >= t_start & time <= t_end;
        event_time = time(event_mask_global);
        event_signal = F_zscore(event_mask_global);
        
        if length(event_signal) < 3
            continue;
        end
        
        % Find local baseline (minimum in window, or start point)
        local_baseline = min(event_signal);
        
        % Calculate AUC from local baseline
        event_signal_corrected = event_signal - local_baseline;
        event_auc = trapz(event_time, event_signal_corrected);
        total_transient_auc = total_transient_auc + event_auc;
        
        % Plot transient window
        event_mask_window = time_window >= t_start & time_window <= t_end;
        if sum(event_mask_window) > 0
            % Shade the transient AUC
            t_plot = time_window(event_mask_window);
            s_plot = signal_window(event_mask_window);
            baseline_plot = ones(size(t_plot)) * local_baseline;
            
            % Fill area between signal and local baseline
            fill([t_plot/60; flipud(t_plot/60)], ...
                 [s_plot; flipud(baseline_plot)], ...
                 colors(i,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
            
            % Mark local baseline
            plot([t_start/60, t_end/60], [local_baseline, local_baseline], ...
                 '--', 'Color', colors(i,:), 'LineWidth', 2);
            
            % Mark peak
            plot(peak_time/60, transients_in_window.peak_amps(i), 'o', ...
                 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), 'MarkerSize', 8);
            
            % Annotate with AUC value
            text(peak_time/60, transients_in_window.peak_amps(i) + 0.5, ...
                 sprintf('AUC=%.1f', event_auc), 'FontSize', 8, 'Color', colors(i,:));
        end
    end
    
    ylabel('Z-score');
    title(sprintf('Event-Based Transient AUC: Each event from its local baseline\nTotal Transient AUC = %.2f', ...
        total_transient_auc), 'FontSize', 11);
    legend({'Baseline (z=0)', 'Signal', 'Transient AUC', 'Local baseline', 'Peak'}, 'Location', 'best');
    grid on;
    xlim([time_window(1)/60, time_window(end)/60]);
    ax = gca;
    ax.XMinorTick = 'on';
    x_start = floor(time_window(1)/60);
    x_end = ceil(time_window(end)/60);
    ax.XAxis.MinorTickValues = x_start:1:x_end;
    
    %% SUBPLOT 3: Sustained AUC (What remains after removing transients)
    subplot(3,1,3);
    hold on;
    
    % Calculate "sustained" component by masking out transient windows
    sustained_signal = signal_window;
    
    for i = 1:length(transients_in_window.peak_times)
        peak_time = transients_in_window.peak_times(i);
        transient_width = transients_in_window.widths(i);
        window_size = max(3.0, transient_width * 2);
        
        t_start = peak_time - window_size;
        t_end = peak_time + window_size;
        
        % Mask out this transient
        mask_out = time_window >= t_start & time_window <= t_end;
        sustained_signal(mask_out) = NaN; % exclude transient periods
    end
    
    % Plot baseline
    plot(time_window/60, zeros(size(time_window)), 'k--', 'LineWidth', 1.5);
    
    % Plot full signal (faded)
    plot(time_window/60, signal_window, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    
    % Plot sustained component (signal with transients removed)
    plot(time_window/60, sustained_signal, 'b-', 'LineWidth', 2);
    
    % Shade positive sustained
    sustained_pos = max(sustained_signal, 0);
    sustained_pos(isnan(sustained_pos)) = 0;
    area(time_window/60, sustained_pos, 'FaceColor', [0.2 0.2 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Calculate sustained AUC
    sustained_auc = trapz(time_window(~isnan(sustained_signal)), ...
                          sustained_signal(~isnan(sustained_signal)));
    
    ylabel('Z-score');
    xlabel('Time (min)');
    title(sprintf('Sustained AUC: Activity after removing transient events\nSustained AUC = %.2f (Net - Transient = %.2f - %.2f)', ...
        sustained_auc, auc_net, total_transient_auc), 'FontSize', 11);
    legend({'Baseline (z=0)', 'Full signal', 'Sustained component', 'Sustained AUC'}, 'Location', 'best');
    grid on;
    xlim([time_window(1)/60, time_window(end)/60]);
    ax = gca;
    ax.XMinorTick = 'on';
    x_start = floor(time_window(1)/60);
    x_end = ceil(time_window(end)/60);
    ax.XAxis.MinorTickValues = x_start:1:x_end;
    
    %% Save figure
    sgtitle(sprintf('AUC Component Diagnostic - %s (%.1f min ± %.1f min)', ...
        mouseName, diagnostic_center_min, window_min), 'FontSize', 14, 'FontWeight', 'bold');
    
    savePath = fullfile(outFolder, sprintf('AUC_diagnostic_%.1fmin_%s.png', diagnostic_center_min, mouseName));
    saveas(fig, savePath);
    close(fig);
    
    fprintf('Saved AUC diagnostic figure: %s\n', savePath);
end

function [total_auc, mean_auc] = calculateEventBasedAUC(F_zscore, time, transients, transient_mask)
    % Calculate AUC for each transient event from its local baseline
    %
    % Inputs:
    %   F_zscore - full z-scored signal
    %   time - full time vector
    %   transients - transient structure with peak_times, locs, widths
    %   transient_mask - logical mask for which transients to include
    %
    % Outputs:
    %   total_auc - sum of all event AUCs
    %   mean_auc - mean AUC per event
    
    if ~any(transient_mask)
        total_auc = 0;
        mean_auc = 0;
        return;
    end
    
    % Get transients in this period
    peak_times_period = transients.peak_times(transient_mask);
    locs_period = transients.locs(transient_mask);
    widths_period = transients.widths(transient_mask);
    
    event_aucs = zeros(length(peak_times_period), 1);
    
    for i = 1:length(peak_times_period)
        peak_time = peak_times_period(i);
        peak_idx = locs_period(i);
        transient_width = widths_period(i);
        
        % Define window around transient (±2 × width, minimum 3 sec on each side)
        window_size = max(3.0, transient_width * 2);
        t_start = peak_time - window_size;
        t_end = peak_time + window_size;
        
        % Extract event window
        event_mask = time >= t_start & time <= t_end;
        
        if sum(event_mask) < 3
            continue; % Skip if too few points
        end
        
        event_time = time(event_mask);
        event_signal = F_zscore(event_mask);
        
        % Find local baseline (minimum in window)
        local_baseline = min(event_signal);
        
        % Calculate AUC from local baseline
        event_signal_corrected = event_signal - local_baseline;
        event_aucs(i) = trapz(event_time, event_signal_corrected);
    end
    
    % Remove any zero entries (failed events)
    event_aucs = event_aucs(event_aucs > 0);
    
    total_auc = sum(event_aucs);
    
    if ~isempty(event_aucs)
        mean_auc = mean(event_aucs);
    else
        mean_auc = 0;
    end
end

function createTimeBinnedSummary(F_zscore, time, transients, timeBins, baseline_stats, outFolder, mouseName, bin_size_min)
    % Create time-binned summary tables post-injection
    %
    % Inputs:
    %   bin_size_min - bin size in minutes (e.g., 60 for hourly, 30 for 30-min)
    
    % Find injection start and end times
    inj_start = timeBins.injStart;
    
    if isempty(timeBins.injEnd)
        warning('No injection end time specified. Using injection start.');
        post_inj_start = inj_start;
    else
        post_inj_start = timeBins.injEnd;
    end
    
    % Define time range for binning (post-injection to end)
    max_time = max(time);
    
    if post_inj_start >= max_time
        warning('Injection end time is beyond data range. Cannot create time-binned summary.');
        return;
    end
    
    % Calculate number of bins
    total_duration_min = (max_time - post_inj_start) / 60;
    n_bins = floor(total_duration_min / bin_size_min);
    
    if n_bins < 1
        warning('Insufficient data for time binning (need at least %d min post-injection)', bin_size_min);
        return;
    end
    
    % Initialize output arrays
    time_labels = cell(n_bins, 1);
    mean_zscore = zeros(n_bins, 1);
    delta_zscore = zeros(n_bins, 1);
    event_count = zeros(n_bins, 1);
    event_cumulative = zeros(n_bins, 1);
    event_freq = zeros(n_bins, 1);
    delta_event_freq = zeros(n_bins, 1);
    
    % Get baseline event frequency for comparison
    baseline_mask = time >= (timeBins.injStart - baseline_stats.duration_min * 60) & time < timeBins.injStart;
    baseline_transient_mask = false(size(transients.peak_times));
    for i = 1:length(transients.peak_times)
        if transients.peak_times(i) >= min(time(baseline_mask)) && transients.peak_times(i) < timeBins.injStart
            baseline_transient_mask(i) = true;
        end
    end
    baseline_event_count = sum(baseline_transient_mask);
    baseline_event_freq = baseline_event_count / baseline_stats.duration_min;
    
    % Process each bin
    cumulative_events = 0;
    
    for i = 1:n_bins
        % Define bin boundaries
        bin_start = post_inj_start + (i-1) * bin_size_min * 60;
        bin_end = post_inj_start + i * bin_size_min * 60;
        
        % Create label
        if bin_size_min >= 60
            hours_start = (i-1) * bin_size_min / 60;
            hours_end = i * bin_size_min / 60;
            time_labels{i} = sprintf('%.1f-%.1f hr', hours_start, hours_end);
        else
            time_labels{i} = sprintf('%d-%d min', (i-1)*bin_size_min, i*bin_size_min);
        end
        
        % Get data in this bin
        bin_mask = time >= bin_start & time < bin_end;
        
        if sum(bin_mask) < 1
            % No data in this bin
            mean_zscore(i) = NaN;
            delta_zscore(i) = NaN;
            event_count(i) = 0;
            event_freq(i) = NaN;
            delta_event_freq(i) = NaN;
            continue;
        end
        
        % Calculate mean z-score
        bin_signal = F_zscore(bin_mask);
        mean_zscore(i) = mean(bin_signal);
        delta_zscore(i) = mean_zscore(i) - 0; % baseline mean is 0 by definition
        
        % Count events in this bin
        events_in_bin = sum(transients.peak_times >= bin_start & transients.peak_times < bin_end);
        event_count(i) = events_in_bin;
        
        % Cumulative events
        cumulative_events = cumulative_events + events_in_bin;
        event_cumulative(i) = cumulative_events;
        
        % Event frequency (events per minute)
        bin_duration_min = bin_size_min;
        event_freq(i) = events_in_bin / bin_duration_min;
        
        % Delta from baseline
        delta_event_freq(i) = event_freq(i) - baseline_event_freq;
    end
    
    % Create output table
    summary_table = table(time_labels, mean_zscore, delta_zscore, ...
        event_count, event_cumulative, event_freq, delta_event_freq, ...
        'VariableNames', {'Time_Window', 'Mean_Zscore', 'Delta_Zscore_From_Baseline', ...
        'Ca_Event_Count', 'Ca_Event_Cumulative', 'Ca_Event_Freq_Per_Min', ...
        'Delta_Event_Freq_From_Baseline'});
    
    % Save to file
    if bin_size_min >= 60
        filename = sprintf('time_binned_hourly_%s.csv', mouseName);
    else
        filename = sprintf('time_binned_%dmin_%s.csv', bin_size_min, mouseName);
    end
    
    writetable(summary_table, fullfile(outFolder, filename));
    fprintf('Saved: %s\n', filename);
    
    % Display summary
    fprintf('\n--- TIME-BINNED SUMMARY (%d-min bins) ---\n', bin_size_min);
    fprintf('Number of bins: %d\n', n_bins);
    fprintf('Baseline event frequency: %.2f events/min\n', baseline_event_freq);
    fprintf('Total events post-injection: %d\n', cumulative_events);
end
