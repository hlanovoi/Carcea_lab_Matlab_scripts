

% This script processes photometry data and performs the following:
% 1. Basic signal processing and z-scoring
% 2. 1-minute binning of baseline, post-injection, and post-DSI periods
% 3. Normalization to both full baseline and last 15-min baseline period
% 4. Saves both raw and binned/normalized data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% YOU MUST FILL THIS OUT FIRST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% directory of the data
theFolder = '/Users/HL801/Desktop/Photometry_data/CeA_RF_4cup_C4_M/10.26.23-LPSTESTDAY/269445/1L_and_1R/';
% name your mouse
mouseName = '1R';
% Is this mouse a sensor mouse (YES) or GCamp (NO)?
SigProcMode = 'NO';

%%%%%%%%%%%%%%
% when to begin the analysis, in seconds
beginHere = 0;
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% When was the mouse injected?
injStart = [2370];
injEnd = [2392];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% When did the direct social interaction take place? 
DSIStart = [];
DSIEnd = [];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% When was cage added? 
Cagestart = [];
CageEnd = [];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% When was cage with food added? 
Cagefoodstart = [];
CagefoodEnd = [];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% When was just food added? 
FoodStart = [];
FoodEnd = [];
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% enter beginning and end time (in seconds) of period(s) you want removed
rmStart = [];
rmEnd = [];
%%%%%%%%%%%%%%

% directory you want data saved
outFolder = [theFolder '/Output_Parsed_By1min__1R__4_10_25/'];
% Create the output directory if it does not exist
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Debug Configuration
DEBUG = true;
SAVE_DIAGNOSTICS = true;
if SAVE_DIAGNOSTICS
    diagnosticsFolder = fullfile(outFolder, 'diagnostics');
    if ~exist(diagnosticsFolder, 'dir')
        mkdir(diagnosticsFolder);
    end
    diary(fullfile(diagnosticsFolder, 'processing_log.txt')); % Log all console output
end


% Define header for CSV files
header = {'Time (s)', 'Time (min)', 'Signal'};

cd(theFolder)

% get the fiber photometry data
fpFile = dir([theFolder '/P*.csv']); % Use a wildcard pattern to match CSV files
if length(fpFile) == 1
    fpFile = fpFile.name;
else
    if isempty(fpFile)
        error('No photometry files found in the specified directory: %s', theFolder);
    else
        error('Multiple photometry files found in the specified directory: %s. Please ensure there is only one CSV file.', theFolder);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the data
fp = readtable(fpFile);
theSig = fp.Region3G(fp.LedState==6);
ref = fp.Region3G(fp.LedState==1);
theTime = fp.Timestamp(fp.LedState==6);
theMin = min([length(theSig), length(ref)]);
theSig = theSig(1:theMin); ref = ref(1:theMin); theTime = theTime(1:theMin);    

%LED states

% 6 is green if using trig 1
% 2 is green if using trig 3
% 4 is red if using trig 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Filter and normalize the time
theTime = theTime - theTime(1);

% truncate the tracking
trackingTime = theTime; 
beginHereTrack = find(trackingTime > beginHere, 1);

% Process the signal
if strcmp(SigProcMode, 'NO') % If GCamp correct to the isosbestic
    temp_x = (1:length(ref));
    fit_exp = fit(temp_x', double(ref), 'exp2');
    fitRef2 = robustfit(fit_exp(temp_x), theSig);
    fitSig = fit_exp(temp_x) * fitRef2(2) + fitRef2(1);
    newSig = theSig ./ fitSig;
else % If sensor mouse, skip isosbestic correction
    newSig = theSig;
end

% Interpolate the signal to match the tracking data
newSig = interp1(linspace(trackingTime(1), max(trackingTime), length(theTime)), newSig, trackingTime);

% Adjust tracking time
trackingTimeAdjusted = trackingTime - beginHere;
trackingTimeAdjusted = trackingTimeAdjusted(beginHereTrack:end);
totalTrackingTimeInSeconds = trackingTimeAdjusted(end);

% Remove bad periods
removeThis = false(size(trackingTime));
for i = 1:length(rmStart)
    removeThis = removeThis | (trackingTime >= rmStart(i) & trackingTime <= rmEnd(i));
end
newSig(removeThis) = [];
trackingTime(removeThis) = [];

% Calculate removed time
removedDurations = rmEnd - rmStart;
totalRemovedTime = sum(removedDurations);
correctedTotalTrackingTimeInSeconds = totalTrackingTimeInSeconds - totalRemovedTime;

% Z-score the photometry signal
newSig = zscore(newSig);

% Downsample the segments 
downsampleFactor = 10;





%% New code to bin the data 

%% Define time periods and their indices
% Find indices for the segments (only if the corresponding times are not empty)
injStartIdx = ~isempty(injStart) * (injStart > 0) * find(trackingTime >= injStart, 1, 'first');
injEndIdx = ~isempty(injEnd) * (injEnd > 0) * find(trackingTime >= injEnd, 1, 'first');

% Define baseline period end (everything before injection starts)
if injStartIdx > 0
    baselineEndIdx = injStartIdx - 1;
else
    baselineEndIdx = length(trackingTime);
end

% Initialize other indices
DSIStartIdx = [];
DSIEndIdx = [];
CagestartIdx = [];
CageEndIdx = [];
CagefoodstartIdx = [];
CagefoodEndIdx = [];
FoodStartIdx = [];
FoodEndIdx = [];

% Find indices for all time periods if they exist
if ~isempty(DSIStart) && DSIStart > 0
    DSIStartIdx = find(trackingTime >= DSIStart, 1, 'first');
end
if ~isempty(DSIEnd) && DSIEnd > 0
    DSIEndIdx = find(trackingTime >= DSIEnd, 1, 'first');
end
if ~isempty(Cagestart) && Cagestart > 0
    CagestartIdx = find(trackingTime >= Cagestart, 1, 'first');
end
if ~isempty(CageEnd) && CageEnd > 0
    CageEndIdx = find(trackingTime >= CageEnd, 1, 'first');
end
if ~isempty(Cagefoodstart) && Cagefoodstart > 0
    CagefoodstartIdx = find(trackingTime >= Cagefoodstart, 1, 'first');
end
if ~isempty(CagefoodEnd) && CagefoodEnd > 0
    CagefoodEndIdx = find(trackingTime >= CagefoodEnd, 1, 'first');
end
if ~isempty(FoodStart) && FoodStart > 0
    FoodStartIdx = find(trackingTime >= FoodStart, 1, 'first');
end
if ~isempty(FoodEnd) && FoodEnd > 0
    FoodEndIdx = find(trackingTime >= FoodEnd, 1, 'first');
else
    FoodEndIdx = length(trackingTime);
end





%% Bin the data

% Define bin size
bin_size_minutes = 1;

%% Process baseline period
fprintf('\nProcessing baseline period...\n');

% First, create the binned data for baseline
[baseline_binned, baseline_bins] = createTimeBins(...
    trackingTime(1:baselineEndIdx), ...
    newSig(1:baselineEndIdx), ...
    bin_size_minutes, ...
    DEBUG);

% Analyze transients for each bin
transient_counts = zeros(1, length(baseline_bins.bin_centers));
transient_aucs = zeros(1, length(baseline_bins.bin_centers));
peak_times_cell = cell(1, length(baseline_bins.bin_centers));

for bin_idx = 1:length(baseline_bins.bin_centers)
    bin_mask = baseline_bins.bin_assignment == bin_idx;
    bin_time = trackingTime(1:baselineEndIdx);
    bin_time = bin_time(bin_mask);
    bin_signal = newSig(1:baselineEndIdx);
    bin_signal = bin_signal(bin_mask);
    
    if ~isempty(bin_signal)
        [count, auc, peaks] = findCalciumTransients(bin_time, bin_signal);
        transient_counts(bin_idx) = count;
        transient_aucs(bin_idx) = auc;
        peak_times_cell{bin_idx} = peaks;
    end
end

% Calculate baseline statistics
baseline_mean = mean(baseline_binned.mean);
baseline_std = mean(baseline_binned.std);
last_baseline_bin_mean = baseline_binned.mean(end);

% Create baseline data structure with proper nesting
baseline_data = struct();
baseline_data.time = trackingTime(1:baselineEndIdx);
baseline_data.signal = newSig(1:baselineEndIdx);
baseline_data.binned = struct();
baseline_data.binned.mean = baseline_binned.mean;
baseline_data.binned.std = baseline_binned.std;
baseline_data.binned.sem = baseline_binned.sem;
baseline_data.binned.n_samples = baseline_binned.n_samples;
baseline_data.binned.bin_info = baseline_bins;
baseline_data.transients = struct(...
    'count', transient_counts, ...
    'auc', transient_aucs, ...
    'peak_times', {peak_times_cell});  % Cell array of peak times per bin

% Save baseline data
saveProcessedData(outFolder, '1_baseline', mouseName, baseline_data, [], true);

baseline_results = baseline_data;  % Store baseline data for reference in later processing




%%

% Process post-injection period if available
if ~isempty(injStartIdx) && ~isempty(injEndIdx) && injStartIdx > 0 && injEndIdx > 0
    % Save injection period data
    injectionData = downsampleAndCombine(trackingTime(injStartIdx:injEndIdx), ...
        newSig(injStartIdx:injEndIdx), downsampleFactor);
    writecell(header, [outFolder '2_injection_' mouseName '.csv']);
    writematrix(injectionData, [outFolder '2_injection_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Process post-injection period
    postInjEndIdx = findNextStartIdx({DSIStartIdx, CagestartIdx, CagefoodstartIdx, FoodStartIdx}, trackingTime);
    
    [postinj_binned, postinj_bins] = createTimeBins(...
        trackingTime(injEndIdx+1:postInjEndIdx), ...
        newSig(injEndIdx+1:postInjEndIdx), ...
        bin_size_minutes, ...
        DEBUG);
    
    % Add transient analysis for each bin
    postinj_transient_counts = zeros(1, length(postinj_bins.bin_centers));
    postinj_transient_aucs = zeros(1, length(postinj_bins.bin_centers));
    postinj_peak_times_cell = cell(1, length(postinj_bins.bin_centers));
    
    for bin_idx = 1:length(postinj_bins.bin_centers)
        bin_mask = postinj_bins.bin_assignment == bin_idx;
        bin_time = trackingTime(injEndIdx+1:postInjEndIdx);
        bin_time = bin_time(bin_mask);
        bin_signal = newSig(injEndIdx+1:postInjEndIdx);
        bin_signal = bin_signal(bin_mask);
        
        if ~isempty(bin_signal)
            [count, auc, peaks] = findCalciumTransients(bin_time, bin_signal);
            postinj_transient_counts(bin_idx) = count;
            postinj_transient_aucs(bin_idx) = auc;
            postinj_peak_times_cell{bin_idx} = peaks;
        end
    end
    
    % Create post-injection data structure
    postinj_data = struct();
    postinj_data.time = trackingTime(injEndIdx+1:postInjEndIdx);
    postinj_data.signal = newSig(injEndIdx+1:postInjEndIdx);
    postinj_data.binned = struct();
    postinj_data.binned.mean = postinj_binned.mean;
    postinj_data.binned.std = postinj_binned.std;
    postinj_data.binned.sem = postinj_binned.sem;
    postinj_data.binned.n_samples = postinj_binned.n_samples;
    postinj_data.binned.bin_info = postinj_bins;
    postinj_data.transients = struct(...
        'count', postinj_transient_counts, ...
        'auc', postinj_transient_aucs, ...
        'peak_times', {postinj_peak_times_cell});
    
    % Save post-injection results with baseline reference
    saveProcessedData(outFolder, '3_post_injection', mouseName, postinj_data, baseline_results, false);
end


% Process DSI period if available
if ~isempty(DSIStartIdx) && ~isempty(DSIEndIdx) && DSIStartIdx > 0 && DSIEndIdx > 0
    % Save DSI period data
    DSIData = downsampleAndCombine(trackingTime(DSIStartIdx:DSIEndIdx), ...
        newSig(DSIStartIdx:DSIEndIdx), downsampleFactor);
    writecell(header, [outFolder '4_DSI_' mouseName '.csv']);
    writematrix(DSIData, [outFolder '4_DSI_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Process post-DSI period
    postDSIEndIdx = findNextStartIdx({CagestartIdx, CagefoodstartIdx, FoodStartIdx}, trackingTime);
    
    % Create binned data for post-DSI period
    [postdsi_binned, postdsi_bins] = createTimeBins(...
        trackingTime(DSIEndIdx+1:postDSIEndIdx), ...
        newSig(DSIEndIdx+1:postDSIEndIdx), ...
        bin_size_minutes, ...
        DEBUG);
    
    % Analyze transients for each post-DSI bin
    postdsi_transient_counts = zeros(1, length(postdsi_bins.bin_centers));
    postdsi_transient_aucs = zeros(1, length(postdsi_bins.bin_centers));
    postdsi_peak_times_cell = cell(1, length(postdsi_bins.bin_centers));
    
    for bin_idx = 1:length(postdsi_bins.bin_centers)
        bin_mask = postdsi_bins.bin_assignment == bin_idx;
        bin_time = trackingTime(DSIEndIdx+1:postDSIEndIdx);
        bin_time = bin_time(bin_mask);
        bin_signal = newSig(DSIEndIdx+1:postDSIEndIdx);
        bin_signal = bin_signal(bin_mask);
        
        if ~isempty(bin_signal)
            [count, auc, peaks] = findCalciumTransients(bin_time, bin_signal);
            postdsi_transient_counts(bin_idx) = count;
            postdsi_transient_aucs(bin_idx) = auc;
            postdsi_peak_times_cell{bin_idx} = peaks;
        end
    end
    
    % Create post-DSI data structure
    postdsi_data = struct();
    postdsi_data.time = trackingTime(DSIEndIdx+1:postDSIEndIdx);
    postdsi_data.signal = newSig(DSIEndIdx+1:postDSIEndIdx);
    postdsi_data.binned = struct();
    postdsi_data.binned.mean = postdsi_binned.mean;
    postdsi_data.binned.std = postdsi_binned.std;
    postdsi_data.binned.sem = postdsi_binned.sem;
    postdsi_data.binned.n_samples = postdsi_binned.n_samples;
    postdsi_data.binned.bin_info = postdsi_bins;
    postdsi_data.transients = struct(...
        'count', postdsi_transient_counts, ...
        'auc', postdsi_transient_aucs, ...
        'peak_times', {postdsi_peak_times_cell});
    
    % Save post-DSI results with baseline reference
    saveProcessedData(outFolder, '5_post_DSI', mouseName, postdsi_data, baseline_results, false);
end

%%

baselineData = downsampleAndCombine(trackingTime(1:baselineEndIdx), newSig(1:baselineEndIdx), downsampleFactor);
writecell(header, [outFolder '1_baseline_' mouseName '.csv']);
writematrix(baselineData, [outFolder '1_baseline_' mouseName '.csv'], 'WriteMode', 'append');

% Injection (if injection times are provided)
if ~isempty(injStartIdx) && ~isempty(injEndIdx) && injStartIdx > 0 && injEndIdx > 0
    injectionData = downsampleAndCombine(trackingTime(injStartIdx:injEndIdx), newSig(injStartIdx:injEndIdx), downsampleFactor);
    writecell(header, [outFolder '2_injection_' mouseName '.csv']);
    writematrix(injectionData, [outFolder '2_injection_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Post-injection
    postInjectionEndIdx = findNextStartIdx({DSIStartIdx, CagestartIdx, CagefoodstartIdx, FoodStartIdx}, trackingTime);
    postInjectionData = downsampleAndCombine(trackingTime(injEndIdx+1:postInjectionEndIdx), newSig(injEndIdx+1:postInjectionEndIdx), downsampleFactor);
    writecell(header, [outFolder '3_post_injection_' mouseName '.csv']);
    writematrix(postInjectionData, [outFolder '3_post_injection_' mouseName '.csv'], 'WriteMode', 'append');

end

% DSI (if DSI times are provided)
if ~isempty(DSIStartIdx) && ~isempty(DSIEndIdx) && DSIStartIdx > 0 && DSIEndIdx > 0
    DSIData = downsampleAndCombine(trackingTime(DSIStartIdx:DSIEndIdx), newSig(DSIStartIdx:DSIEndIdx), downsampleFactor);
    writecell(header, [outFolder '4_DSI_' mouseName '.csv']);
    writematrix(DSIData, [outFolder '4_DSI_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Post-DSI
    postDSIEndIdx = findNextStartIdx({CagestartIdx, CagefoodstartIdx, FoodStartIdx}, trackingTime);
    postDSIData = downsampleAndCombine(trackingTime(DSIEndIdx+1:postDSIEndIdx), newSig(DSIEndIdx+1:postDSIEndIdx), downsampleFactor);
    writecell(header, [outFolder '5_post_DSI_' mouseName '.csv']);
    writematrix(postDSIData, [outFolder '5_post_DSI_' mouseName '.csv'], 'WriteMode', 'append');

end


% Cage (if Cage times are provided)
if ~isempty(CagestartIdx) && ~isempty(CageEndIdx) && CagestartIdx > 0 && CageEndIdx > 0
    CageData = downsampleAndCombine(trackingTime(CagestartIdx:CageEndIdx), newSig(CagestartIdx:CageEndIdx), downsampleFactor);
    writecell(header, [outFolder '6_Cage_' mouseName '.csv']);
    writematrix(CageData, [outFolder '6_Cage_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Post-Cage
    postCageEndIdx = findNextStartIdx({CagefoodstartIdx, FoodStartIdx}, trackingTime);
    postCageData = downsampleAndCombine(trackingTime(CageEndIdx+1:postCageEndIdx), newSig(CageEndIdx+1:postCageEndIdx), downsampleFactor);
    writecell(header, [outFolder '7_post_Cage_' mouseName '.csv']);
    writematrix(postCageData, [outFolder '7_post_Cage_' mouseName '.csv'], 'WriteMode', 'append');
end

% Cagefood (if Cagefood times are provided)
if ~isempty(CagefoodstartIdx) && ~isempty(CagefoodEndIdx) && CagefoodstartIdx > 0 && CagefoodEndIdx > 0
    CagefoodData = downsampleAndCombine(trackingTime(CagefoodstartIdx:CagefoodEndIdx), newSig(CagefoodstartIdx:CagefoodEndIdx), downsampleFactor);
    writecell(header, [outFolder '8_Cagefood_' mouseName '.csv']);
    writematrix(CagefoodData, [outFolder '8_Cagefood_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Post-Cagefood
    postCagefoodEndIdx = findNextStartIdx({FoodStartIdx}, trackingTime);
    postCagefoodData = downsampleAndCombine(trackingTime(CagefoodEndIdx+1:postCagefoodEndIdx), newSig(CagefoodEndIdx+1:postCagefoodEndIdx), downsampleFactor);
    writecell(header, [outFolder '9_post_Cagefood_' mouseName '.csv']);
    writematrix(postCagefoodData, [outFolder '9_post_Cagefood_' mouseName '.csv'], 'WriteMode', 'append');
end

% Food (if Food start time is provided)
if ~isempty(FoodStartIdx) && FoodStartIdx > 0
    FoodData = downsampleAndCombine(trackingTime(FoodStartIdx:FoodEndIdx), newSig(FoodStartIdx:FoodEndIdx), downsampleFactor);
    writecell(header, [outFolder '10_Food_' mouseName '.csv']);
    writematrix(FoodData, [outFolder '10_Food_' mouseName '.csv'], 'WriteMode', 'append');
    
    % Post-Food is not needed in this case as Food extends to the end of recording
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [downsampledData] = downsampleAndCombine(time, signal, factor)
    downsampledTime = downsample(time, factor);
    downsampledSignal = downsample(signal, factor);
    downsampledMinutes = downsampledTime / 60; % Convert seconds to minutes
    downsampledData = [downsampledTime, downsampledMinutes, downsampledSignal];
end

function nextIdx = findNextStartIdx(indices, trackingTime)
    nextIdx = find(~cellfun(@isempty, indices), 1);
    if isempty(nextIdx)
        nextIdx = length(trackingTime);
    else
        nextIdx = indices{nextIdx}(1) - 1;
    end
end


function normalized_data = normalizeToReference(data, reference_value)
    normalized_data = data ./ reference_value - 1;  % Convert to % change from baseline
end

function verifyBaselinePeriod(time, injStart)
    if ~isempty(injStart) && injStart > 0
        actual_duration = time(end) - time(1);
        expected_duration = injStart - time(1);
        if abs(actual_duration - expected_duration) > 1  % 1 second tolerance
            warning('Baseline period duration mismatch. Check injection start time.');
        end
    end
end

function [num_transients, total_auc, peak_times] = findCalciumTransients(time, signal)
    % Ensure time and signal are column vectors
    time = time(:);
    signal = signal(:);
    
    % Parameters optimized for your data
    min_peak_prominence = 0.5;     % Set to detect peaks above 0.5 z-scores
    min_peak_distance = 0.2;       % 200ms minimum between peaks
    
    % Simplified kinetics parameters
    min_rise_decay_ratio = 1.2;    % Modest requirement for rise vs decay
    max_rise_duration = 0.3;       % Allow slightly longer rise times
    min_decay_duration = 0.2;      % Allow shorter decay phases
    
    % Initialize output variables
    num_transients = 0;
    total_auc = 0;
    peak_times = [];
    
    % Check if we have enough data points
    if length(signal) < 3
        return;
    end
    
    % Find peaks with simplified criteria
    [peaks, peak_locs, widths] = findpeaks(signal, ...
        'MinPeakProminence', min_peak_prominence * std(signal), ...
        'MinPeakDistance', round(min_peak_distance * (length(time)/max(time))));
    
    % Process each detected peak
    for j = 1:length(peaks)
        peak_idx = peak_locs(j);
        
        % Define analysis windows
        half_width = round(widths(j)/2);
        start_idx = max(1, peak_idx - half_width);
        end_idx = min(length(signal), peak_idx + half_width);
        
        % Calculate rise and decay windows
        rise_window = start_idx:peak_idx;
        decay_window = peak_idx:end_idx;
        
        if length(rise_window) >= 2 && length(decay_window) >= 2
            % Basic kinetics check
            rise_slope = diff(signal(rise_window))./diff(time(rise_window));
            decay_slope = diff(signal(decay_window))./diff(time(decay_window));
            
            mean_rise_slope = mean(rise_slope);
            mean_decay_slope = mean(decay_slope);
            
            if mean_rise_slope > 0 && mean_decay_slope < 0
                % Valid transient detected
                num_transients = num_transients + 1;
                peak_times = [peak_times; time(peak_idx)];
                
                % Calculate AUC
                baseline = min(signal(start_idx:end_idx));
                transient_signal = signal(start_idx:end_idx) - baseline;
                transient_time = time(start_idx:end_idx);
                total_auc = total_auc + trapz(transient_time, transient_signal);
            end
        end
    end
end





%%
% NEW FUNCTIONS
%%

function normalized_data = normalizeToBaseline(data, baseline_value)
    % Normalize data relative to baseline using (Value - Baseline)/Baseline
    normalized_data = (data - baseline_value) ./ baseline_value;
end

function [diff_counts, diff_auc] = calculateTransientDifferences(current_counts, current_auc, baseline_counts, baseline_auc)
    % Calculate differences in transient counts and AUC from baseline
    diff_counts = current_counts - baseline_counts;
    diff_auc = current_auc - baseline_auc;
end

function [transient_data] = analyzeTransients(time, signal, bin_info, DEBUG)
    transient_data = struct();
    transient_data.count = zeros(1, length(bin_info.bin_centers));
    transient_data.auc = zeros(1, length(bin_info.bin_centers));
    
    for i = 1:length(bin_info.bin_centers)
        bin_mask = bin_info.bin_assignment == i;
        bin_time = time(bin_mask);
        bin_signal = signal(bin_mask);
        
        if ~isempty(bin_signal)
            [num_trans, auc] = findCalciumTransients(bin_time, bin_signal);
            transient_data.count(i) = num_trans;
            transient_data.auc(i) = auc;
            
            if DEBUG
                fprintf('Bin %d: %d transients, AUC = %.2f\n', ...
                    i, num_trans, auc);
            end
        end
    end
    
    if DEBUG
        fprintf('\nTransient Analysis Summary:\n');
        fprintf('Total transients detected: %d\n', sum(transient_data.count));
        fprintf('Mean transients per bin: %.2f\n', mean(transient_data.count));
    end
end


function saveProcessedData(outFolder, prefix, mouseName, data, reference_data, is_baseline)
    % Create filename for raw data
    raw_filename = fullfile(outFolder, sprintf('%s_raw_%s.csv', prefix, mouseName));
    
    % Create raw data table
    raw_table = table(data.time, data.signal, 'VariableNames', {'Time_sec', 'Signal'});
    
    % Add bin assignment if available
    if isfield(data.binned, 'bin_info') && isfield(data.binned.bin_info, 'bin_assignment')
        raw_table.Bin_Number = data.binned.bin_info.bin_assignment;
    end
    
    % Save raw data
    writetable(raw_table, raw_filename);
    
    % Create binned data filename
    binned_filename = fullfile(outFolder, sprintf('%s_binned_%s.csv', prefix, mouseName));
    
    % Process binned data if available
    if isfield(data, 'binned') && isfield(data.binned, 'bin_info')
         % Convert bin centers to minutes and adjust labels to show end of bin time
        bin_centers = data.binned.bin_info.bin_centers;
        bin_size_minutes = 1;  % Using the known bin size
        relative_time = ((bin_centers - bin_centers(1))/60) + bin_size_minutes;  % Add bin_size to show end time
        num_bins = length(bin_centers);
        
        % Initialize the table with the basic binned data
        binned_table = table(...
            bin_centers', ... % Time in seconds
            relative_time', ... % Time in minutes relative to start
            data.binned.mean', ... % Mean signal
            data.binned.std', ... % Standard deviation
            data.binned.sem', ... % Standard error of mean
            data.binned.n_samples', ... % Number of samples
            'VariableNames', {...
                'Time_sec', ...
                'Time_relative_min', ...
                'Signal_mean', ...
                'Signal_std', ...
                'Signal_sem', ...
                'Samples_per_bin'});
        
        % Add delta F/F calculations if reference data is provided
        if ~is_baseline && ~isempty(reference_data)
            % Calculate baseline means across all baseline bins
            baseline_mean_signal = mean(reference_data.binned.mean);
            
            % Calculate delta F/F for mean signal
            delta_f_f_mean = (data.binned.mean - baseline_mean_signal) ./ baseline_mean_signal;
            binned_table.DeltaF_F_mean = delta_f_f_mean';
            
            % Calculate delta F/F for AUC if transient data is available
            if isfield(data, 'transients') && isfield(data.transients, 'auc')
                baseline_mean_auc = mean(reference_data.transients.auc);
                if baseline_mean_auc ~= 0  % Prevent division by zero
                    delta_f_f_auc = (data.transients.auc - baseline_mean_auc) ./ baseline_mean_auc;
                    binned_table.DeltaF_F_AUC = delta_f_f_auc';
                end
            end
        end
        
        % Add transient data if available
        if isfield(data, 'transients')
            if isfield(data.transients, 'count')
                binned_table.Transient_count = data.transients.count';
            end
            if isfield(data.transients, 'auc')
                binned_table.Transient_AUC = data.transients.auc';
            end
        end
        
        % Save binned data
        writetable(binned_table, binned_filename);
    end
end

%% Modified createTimeBins function
function [binned_data, bin_info] = createTimeBins(time, signal, bin_size_minutes, DEBUG)
    % Calculate bin parameters
    bin_size_seconds = bin_size_minutes * 60;
    start_time = floor(min(time) / bin_size_seconds) * bin_size_seconds;
    end_time = ceil(max(time) / bin_size_seconds) * bin_size_seconds;
    
    bin_edges = start_time:bin_size_seconds:end_time;
    num_bins = length(bin_edges) - 1;
    
    % Initialize output structures
    binned_data = struct();
    binned_data.mean = zeros(1, num_bins);
    binned_data.std = zeros(1, num_bins);
    binned_data.sem = zeros(1, num_bins);
    binned_data.n_samples = zeros(1, num_bins);
    
    bin_info = struct();
    bin_info.bin_edges = bin_edges;
    bin_info.bin_centers = bin_edges(1:end-1) + bin_size_seconds/2;
    bin_info.bin_assignment = zeros(size(time));
    
    % Process each bin
    for i = 1:num_bins
        bin_mask = time >= bin_edges(i) & time < bin_edges(i+1);
        bin_data = signal(bin_mask);
        bin_info.bin_assignment(bin_mask) = i;
        
        if ~isempty(bin_data)
            binned_data.mean(i) = mean(bin_data);
            binned_data.std(i) = std(bin_data);
            binned_data.sem(i) = std(bin_data)/sqrt(length(bin_data));
            binned_data.n_samples(i) = length(bin_data);
        end
    end
    
    if DEBUG
        fprintf('Binning Analysis:\n');
        fprintf('Total time range: %.2f to %.2f seconds\n', min(time), max(time));
        fprintf('Number of bins: %d\n', num_bins);
        fprintf('Mean samples per bin: %.2f\n', mean(binned_data.n_samples));
    end
end
