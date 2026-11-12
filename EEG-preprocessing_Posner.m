%-----------------------------------------------------------------------%
% EEG preprocessing pipeline for the Posner spatial attention task
%-----------------------------------------------------------------------%
% uses: MATLAB R2025a | EEGLAB 2025.0.0

% Author: Leonie Polzer

% EEG data structure:
% - Participants: children/adolescents
% - Recording: 64 EEG channels (FCz reference) + 1 vertical EOG channel;
%   sampling rate 1000 Hz
% - Task: Posner spatial cueing task with spatial cues (left/right, 50 ms) followed by targets
%   (left/right, 500 ms) after 50 or 200 ms
% - Conditions: Valid cue + target, invalid cue + target,
%   cue only (catch), or target only (uncued)
% - Trials: 640 total (356 valid, 156 invalid, 64 catch, 64 uncued).
% - Responses: Left/right keypress according to target location
% - EEG events: Trigger codes distinguish cue/target location,
%   validity, stimulus timing, and behavioral responses

% Description:
% - Loads and preprocesses continuous EEG data, including filtering, bad
%   channel detection, re-referencing, and ICA-based artifact correction.
% - Creates epochs by condition and stimulus location,
%   performs baseline correction and semi-automated trial rejection,
%   and prepares cleaned data for time-frequency analyses
% - Processes behavioral responses

%% Install session
addpath('...\eeglab2025.0.0')

% open eeglab
eeglab

% keep figures from docking/stacking
set(0,'DefaultFigureWindowStyle','normal'); 

% data path
% data_folder = 'S:\KJP_Neurophys\Studien\STIPED\WP5 Querschnitt (SoKoASS)\6_Projekte\1_Posner\2_Paper Leonie - oscillations\test_dataset_2\';
data_folder = '...\Posner\raw files\'

file_list = dir(fullfile(data_folder, '*.vhdr'));

% empty array for all eeg-datasets
EEG = [];

for i = 1:length(file_list)

    % load BrainVision-file
    EEG = pop_loadbv(data_folder, file_list(i).name);
    
    % file name without ending
    [~, name, ~] = fileparts(file_list(i).name);
    
    % set file name
    EEG.setname = name;
    
    % take file name as info
    EEG.filename = file_list(i).name;

    % Extract subject name
    idxT2     = strfind(upper(name), '_T2');
    idxPosner = strfind(upper(name), '_POSNER');
    
    if ~isempty(idxT2)
        EEG.subject = name(1:idxT2(1)-1);
    elseif ~isempty(idxPosner)
        EEG.subject = name(1:idxPosner(1)-1);
    else
        EEG.subject = ['subj' num2str(i)];  % fallback
    end

    % assign group according to file
    upper_name = upper(name);  % upper case

    if contains(upper_name, 'ASS') || contains(upper_name, 'AUT') || contains(upper_name, 'QS') % remove QS if not SOKO study
        EEG.group = 'ASD';
    elseif contains(upper_name, 'HC')
        EEG.group = 'NT';
    else
        EEG.group = 'Unknown';  % fallback
    end

    % ---------------------------
    % 1. Load Channel Locations
    % ---------------------------
    chanlocs_file = '...\64_extended20-10.elp';
    
    % % Register FCz as recording reference
    EEG = pop_chanedit(EEG, ...
        'load', {chanlocs_file, 'filetype', 'autodetect'}, ...
        'changefield', {65, 'datachan', 0}, ...
        'setref', {'1:65', 'FCz'});
    EEG = eeg_checkset(EEG);

    % ---------------------------
    % 2. Downsampling to 500 Hz
    % ---------------------------
    EEG = pop_resample(EEG, 500);
    EEG = eeg_checkset(EEG);

    % ---------------------------
    % 3. Remove line noise 
    % ---------------------------
    EEG        = pop_cleanline(EEG, ...
    'bandwidth',        2, ...
    'chanlist',         1:64, ...
    'computepower',     1, ...
    'linefreqs',        [50 100], ...
    'newversion',       0, ...
    'normSpectrum',     0, ...
    'p',                0.01, ...
    'pad',              2, ...
    'plotfigures',      0, ...
    'scanforlines',     1, ...
    'sigtype',          'Channels', ...
    'taperbandwidth',   2, ...
    'tau',              100, ...
    'verb',             1, ...
    'winsize',          4, ...
    'winstep',          2);

    % -------------------------------
    % 4. Identify bad channels
    % -------------------------------
    
    % 5.1. Before CleanRawData: EEG-Backup
    EEG_orig = EEG;  % copy data set with all channels
    all_labels_orig = {EEG_orig.chanlocs.labels};
    
    % 5.2. CleanRawData (ignore IO)
    % (deletes bad channels, but only in EEG_clean data set)
    EEG_clean = pop_clean_rawdata(EEG, ...
        'FlatlineCriterion', 5, ...
        'ChannelCriterion', 0.65, ...
        'LineNoiseCriterion', 5, ...
        'Highpass', 'off', ...
        'BurstCriterion', 'off', ...
        'WindowCriterion', 'off', ...
        'BurstRejection', 'off', ...
        'Distance', 'Euclidian', ...
        'channels_ignore', {'IO'});
    
    EEG_clean = eeg_checkset(EEG_clean);
    
    % 5.3. Identify excluded channels 
    % channels of channel-cleaned data set (does not include IO, see above)
    labels_clean = {EEG_clean.chanlocs.labels};
    % all original channels (excluding IO)
    labels_orig_noIO = setdiff(all_labels_orig, {'IO'});
    % compare data sets to print channels that would be excluded by algorithm
    bad_chans = setdiff(labels_orig_noIO, labels_clean);
    
    fprintf('\n Bad channels as identified by algorithm:\n');
    disp(bad_chans');

    % 5.4. Visual inspection of channels
    fprintf(' Open EEG-Browser for visual inspection...\n');
    % --- Continuous scroll viewer (own window)
    hScroll = eegplot(EEG_orig.data, ...
        'srate', EEG_orig.srate, ...
        'winlength', 8, ...                              % seconds per screen
        'dispchans', size(EEG_orig.data,1), ...
        'spacing', 50, ...
        'eloc_file', EEG_orig.chanlocs, ...
        'events', EEG_orig.event, ...
        'title', sprintf('%s | Continuous data (original)', EEG_orig.setname));
    % optional: turn on DC offset removal in GUI for better inspection: Mind that this helps
    % identify channel trends but may cover deflections
    
    % --- Spectral overview in a different figure
    plotchans = find(~strcmpi({EEG_orig.chanlocs.labels}, 'IO'));   % ignore IO
    hSpect = figure('Name', sprintf('%s | Spectral overview', EEG_orig.setname), ...
                    'NumberTitle','off');
    clf(hSpect);  % ensure clean canvas
    pop_spectopo(EEG_orig, 1, [EEG_orig.xmin*1000 EEG_orig.xmax*1000], 'EEG', ...
        'freqrange',[2 70], 'electrodes','off', 'plotchans', plotchans);

    % add for spectral overview cheat sheet translating channel numbers
    % into channel labels
    allLabels = {EEG_orig.chanlocs.labels};
    badSet = [];
    if exist('bad_chans','var') && ~isempty(bad_chans)
        badSet = string(bad_chans);   % auto-flagged (from CleanRawData)
    end
    
    % Build mapping "index : label  (flag)" for the plotted channels
    mapStr = strings(numel(plotchans),1);
    for k = 1:numel(plotchans)
        idx = plotchans(k);
        lab = allLabels{idx};
        flag = "";
        if ~isempty(badSet) && any(strcmpi(badSet, lab))
            flag = "  ⟵ auto-flag";
        end
        mapStr(k) = sprintf('%3d : %s%s', idx, lab, flag);
    end
    
    % Put a scrollable listbox on the right side of the spectopo figure
    uicontrol('Parent', hSpect, 'Style','listbox', ...
        'Units','normalized', 'Position',[0.82 0.10 0.17 0.80], ...
        'String', cellstr(mapStr), ...
        'FontName','Consolas', 'FontSize', 10, ...
        'Tag','chanmap');
    
    % Nudge existing axes left to make space for the listbox 
    axs = findall(hSpect, 'Type','axes');
    for ax = transpose(axs)
        set(ax,'Units','normalized');
        pos = get(ax,'Position');
        pos(3) = pos(3)*0.80;  % shrink width by 20%
        set(ax,'Position',pos);
    end
    
    % Place windows side-by-side on Windows
    try
        movegui(hScroll, 'west'); 
        movegui(hSpect,  'east');
    end
    
    % Break for visual inspection
    uiwait(msgbox('Please review electrode activity in both windows. Click OK, if you want to continue.', 'Manual Inspection', 'modal'));


    % ---------------------------
    % 5.5 Manual channel selection (GUI) + store audit
    % ---------------------------
    all_labels_orig = {EEG_orig.chanlocs.labels};
    
    % Open channel picker GUI; manual selection of channels to remove
    [selLabels, selIdx] = pop_chansel(all_labels_orig, ...
                                      'withindex','on', ...
                                      'selectionmode','multiple');
    
    % --- Normalize GUI outputs
    if ischar(selLabels)                      
        selLabels = {selLabels};
    end
    if isempty(selLabels) && ~isempty(selIdx) 
        if islogical(selIdx), selIdx = find(selIdx); end
        if iscell(selIdx) && all(cellfun(@(c) ~isempty(regexp(c,'^\d+$','once')), selIdx))
            selIdx = cellfun(@str2double, selIdx);
        end
        if isnumeric(selIdx) && ~isempty(selIdx)
            selLabels = all_labels_orig(selIdx);
        end
    end
    if isstring(selLabels), selLabels = cellstr(selLabels); end
    
    % --- Apply removal (if any)
    if ~isempty(selLabels)
        EEG = pop_select(EEG_orig, 'nochannel', selLabels);   % remove the selected labels
        rem_chans = setdiff(all_labels_orig, {EEG.chanlocs.labels});
    else
        EEG = EEG_orig;                                       % user cancelled or selected none
        rem_chans = {};
        fprintf('No channels selected for removal.\n');
    end
    
    EEG = eeg_checkset(EEG);
    
    % --- Store bad-channel audit (removed + suggested)
    if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc), EEG.etc = struct; end
    
    % 'bad_chans' from yCleanRawData suggestion earlier:
    if exist('bad_chans','var')
        if isstring(bad_chans), bad_chans = cellstr(bad_chans); end
        if ischar(bad_chans),  bad_chans  = {bad_chans};       end
        % keep only labels that actually exist in the original montage
        bad_chans = intersect(bad_chans, all_labels_orig, 'stable');
    else
        bad_chans = {};
    end
    
    EEG.etc.badchan_suggested = bad_chans;    % from auto algorithm (for traceability)
    EEG.etc.badchan_labels    = rem_chans;    % actually removed (manual)
    EEG.etc.badchan_method    = 'manual_via_pop_chansel';
    EEG.etc.badchan_timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

  
    % ---------------------------
    % 5. Re-referencing to average
    % ---------------------------

    % channels that should not be calculated into average reference
    channels_to_exclude = {'IO'};  
    
    % get indices of these channels
    exclude_inds = find(ismember({EEG.chanlocs.labels}, channels_to_exclude));
    
    % Average-Re-referencing without interpolation or reference retrieval
    disp('>> before reref');
    EEG = pop_reref(EEG, [], ...
        'exclude', exclude_inds, ...
        'interpchan', 'off'); 
    disp('>> after reref');

    EEG = eeg_checkset(EEG);


    % ---------------------------
    % 6. High-pass filter
    % ---------------------------
    % create temporary data set used for ICA
    EEG_for_ICA      = EEG;     
    
    % 0.1 Hz high-pass filter for data
    EEG = pop_eegfiltnew(EEG, 0.1, []);
    EEG = eeg_checkset(EEG);

    % 1 Hz high-pass filter for ICA
    EEG_for_ICA = pop_eegfiltnew(EEG_for_ICA, 1, []);
    EEG_for_ICA = eeg_checkset(EEG_for_ICA);


    % ---------------------------
    % 7. Artifact removal
    % ---------------------------
    % --- 1) Auto-detect on COPY; real data stays untouched
    idxKeep = find(~strcmpi({EEG_for_ICA.chanlocs.labels}, 'IO'));  % exclude IO from detection
    EEG_tmp = EEG_for_ICA;
    
    [~, winrej_auto] = pop_rejcont(EEG_tmp, ...
        'elecrange', idxKeep, ...
        'threshold', 5, ...           
        'freqlimit', [30 80], ...    
        'epochlength', 1, 'overlap', 0.5, ...
        'contiguous', 2, 'addlength', 0.25, 'taper', 'hamming', ...
        'eegplot', 'off', ...
        'onlyreturnselection', 'on', ...
        'verbose', 'off');
    
    % --- 2) Open the editor with auto-marks preloaded 
    nch = size(EEG_for_ICA.data, 1);                 % show all channels
    plotTitle = sprintf('%s | Edit marked windows for ICA', EEG_for_ICA.setname);
    
    eegplot(EEG_for_ICA.data, ...
        'srate', EEG_for_ICA.srate, ...
        'winlength', 8, ...
        'dispchans', nch, ...
        'spacing', 50, ...                            
        'eloc_file', EEG_for_ICA.chanlocs, ...
        'events', EEG_for_ICA.event, ...
        'winrej', winrej_auto, ...
        'butlabel', 'Save selection & close', ...
        'title', plotTitle, ...
        'command', 'assignin(''base'',''EEG_WINREJ_FINAL'', TMPREJ); close(gcbf);');
    
    % --- 3) Wait until the editor window closes (user clicked the button)
    drawnow; pause(0.1);
    hFig = findall(0,'type','figure','tag','EEGPLOT','Name',plotTitle);
    if isempty(hFig)
        figs = findall(0,'type','figure','tag','EEGPLOT');
        if ~isempty(figs), hFig = figs(1); end
    end
    if ~isempty(hFig) && ishghandle(hFig)
        waitfor(hFig);   % returns when the editor is closed
    end
    
    % --- 4) Retrieve edited selection
    winrej_final = [];
    if evalin('base','exist(''EEG_WINREJ_FINAL'',''var'')')
        winrej_final = evalin('base','EEG_WINREJ_FINAL');
        evalin('base','clear EEG_WINREJ_FINAL');
    end
    if ~isempty(winrej_final) && size(winrej_final,2) > 2
        winrej_final = winrej_final(:,1:2);  % keep [start end] sample columns
    end
    
    % --- 5) Optional: confirm and actually cut (APPLIES TO ICA DATASET)
    if ~isempty(winrej_final)
        total_sec = sum((winrej_final(:,2)-winrej_final(:,1)+1))/EEG_for_ICA.srate;
        reply = questdlg(sprintf('Cut %d windows (%.1f s total) from ICA dataset?', ...
                                 size(winrej_final,1), total_sec), ...
                         'Confirm cut', 'Yes','No','Yes');
        if strcmpi(reply,'Yes')
            precut_pnts = EEG_for_ICA.pnts;
            EEG_for_ICA = eeg_eegrej(EEG_for_ICA, winrej_final);   % <-- cut ONCE
            EEG_for_ICA = eeg_checkset(EEG_for_ICA);
    
            % --- audit summary after cutting ---
            if ~isfield(EEG_for_ICA,'etc') || ~isstruct(EEG_for_ICA.etc)
                EEG_for_ICA.etc = struct;
            end
            EEG_for_ICA.etc.winrej_applied = winrej_final;   % keep what you applied
    
            samples_removed_calc = sum(winrej_final(:,2) - winrej_final(:,1) + 1);
            sec_removed_calc     = samples_removed_calc / EEG_for_ICA.srate;
            delta_pnts           = precut_pnts - EEG_for_ICA.pnts;
            n_boundaries         = sum(strcmpi({EEG_for_ICA.event.type}, 'boundary'));
    
            fprintf('\n[ICA cut summary]\n');
            fprintf('  Windows removed:   %d\n', size(winrej_final,1));
            fprintf('  Seconds removed:   %.2f s (planned)\n', sec_removed_calc);
            fprintf('  Samples removed:   %d (planned)\n', samples_removed_calc);
            fprintf('  Samples removed:   %d (delta in data)\n', delta_pnts);
            fprintf('  Boundary events:   %d\n', n_boundaries);
            if abs(delta_pnts - samples_removed_calc) <= 2
                fprintf('  Check: OK (delta ~= planned)\n');
            else
                fprintf('  Check: MISMATCH (investigate rounding/overlap)\n');
            end
        end
    end
    
    % % --- 6) Visual check after rejection (shows boundary events)
    eegplot(EEG_for_ICA.data, ...
        'srate', EEG_for_ICA.srate, ...
        'winlength', 8, ...
        'dispchans', size(EEG_for_ICA.data,1), ...
        'spacing', 50, ...
        'eloc_file', EEG_for_ICA.chanlocs, ...
        'events', EEG_for_ICA.event, ...
        'title', sprintf('%s | Post-cut (boundary events shown)', EEG_for_ICA.setname));
    
    % After applying/removing channels:
    EEG = eeg_checkset(EEG);
    
    % Robust close of all viewers
    try
        % Any eegplot windows 
        hEEG = findall(0,'type','figure','tag','EEGPLOT');
        if ~isempty(hEEG), close(hEEG); end
    
        % Our spectral overview window (by exact name)
        hSpecList = findall(0,'type','figure','Name', sprintf('%s | Spectral overview', EEG_orig.setname));
        if ~isempty(hSpecList), close(hSpecList); end
    end

    
    % --- 7) Write audit log to TXT
    out_folder = '...\Pre_ICA_artifactdetection';
    if ~exist(out_folder, 'dir'), mkdir(out_folder); end
    
    if ~exist('winrej_final','var') || isempty(winrej_final)
        n_windows = 0; samples_removed_planned = 0; sec_removed_calc = 0;
    else
        n_windows = size(winrej_final,1);
        samples_removed_planned = sum(winrej_final(:,2) - winrej_final(:,1) + 1);
        sec_removed_calc = samples_removed_planned / EEG_for_ICA.srate;
    end
    if exist('precut_pnts','var')
        samples_removed_delta = max(precut_pnts - EEG_for_ICA.pnts, 0);
    else
        samples_removed_delta = 0;
    end
    n_boundaries = sum(strcmpi({EEG_for_ICA.event.type}, 'boundary'));
    
    % Build the struct (and keep a copy inside EEG_for_ICA)
    if ~isfield(EEG_for_ICA, 'etc') || ~isstruct(EEG_for_ICA.etc), EEG_for_ICA.etc = struct; end
    participantVal = EEG_for_ICA.subject;
    if isempty(participantVal), participantVal = EEG_for_ICA.setname; end
    
    EEG_for_ICA.etc.artifact_cut_log = struct( ...
        'participant',               participantVal, ...
        'dataset',                   EEG_for_ICA.setname, ...
        'srate',                     EEG_for_ICA.srate, ...
        'n_windows',                 n_windows, ...
        'sec_removed',               sec_removed_calc, ...
        'samples_removed_planned',   samples_removed_planned, ...
        'samples_removed_delta',     samples_removed_delta, ...
        'n_boundaries',              n_boundaries, ...
        'timestamp',                 datestr(now, 'yyyy-mm-dd HH:MM:SS') ...
    );
    
    % Make a safe filename from participant
    subj = EEG_for_ICA.etc.artifact_cut_log.participant;
    if isempty(subj), subj = EEG_for_ICA.setname; end
    safeSubj = regexprep(subj, '[<>:"/\\|?*]+', '_');
    safeSubj = regexprep(safeSubj, '[\. ]+$', '');
    out_file = fullfile(out_folder, [safeSubj '_pre_ICA_artifactdetection.txt']);
    
    % Write as simple key: value text
    fid = fopen(out_file, 'w'); assert(fid ~= -1, 'Could not open file: %s', out_file);
    fprintf(fid, 'participant: %s\n', EEG_for_ICA.etc.artifact_cut_log.participant);
    fprintf(fid, 'dataset: %s\n',     EEG_for_ICA.etc.artifact_cut_log.dataset);
    fprintf(fid, 'srate_Hz: %.6g\n',  EEG_for_ICA.etc.artifact_cut_log.srate);
    fprintf(fid, 'n_windows: %d\n',   EEG_for_ICA.etc.artifact_cut_log.n_windows);
    fprintf(fid, 'sec_removed: %.3f\n', EEG_for_ICA.etc.artifact_cut_log.sec_removed);
    fprintf(fid, 'samples_removed_planned: %d\n', EEG_for_ICA.etc.artifact_cut_log.samples_removed_planned);
    fprintf(fid, 'samples_removed_delta: %d\n',   EEG_for_ICA.etc.artifact_cut_log.samples_removed_delta);
    fprintf(fid, 'n_boundaries: %d\n', EEG_for_ICA.etc.artifact_cut_log.n_boundaries);
    fprintf(fid, 'timestamp: %s\n',    EEG_for_ICA.etc.artifact_cut_log.timestamp);
    fclose(fid);
    fprintf('[Log] Wrote artifact-cut summary to:\n  %s\n', out_file);

    % ---------------------------
    % 8. Save data sets and end loop
    % ---------------------------
    % allows running all the ICAs together in the next step
    EEG_for_ICA = eeg_checkset(EEG_for_ICA);
    
    % Save both datasets with clear names
    outA = '...\Pre_ICA_Sets';
    if ~exist(outA,'dir'), mkdir(outA); end
    
    subj = EEG.subject; if isempty(subj), subj = EEG.setname; end
    base = regexprep(subj,'[<>:"/\\|?*]+','_');
    
    pop_saveset(EEG_for_ICA, 'filename', [base '_preICA_forICA.set'],     'filepath', outA);
    pop_saveset(EEG,          'filename', [base '_preICA_analysis.set'],  'filepath', outA);

    % right before 'end' of the for-loop:
    clear EEG_clean EEG_orig EEG_tmp labels_clean bad_chans rem_chans winrej_auto winrej_final hScroll hSpect
end




%%
 % % ---------------------------
 % % 9. Run ICA 
 % % ---------------------------
 % ICA on all pre-ICA sets

inA  = '...\Pre_ICA_Sets';
outB = '...\ICA_Results';

if ~isfolder(outB), mkdir(outB); end

files = dir(fullfile(inA, '*_preICA_forICA.set'));
fprintf('[Batch ICA] Found %d pre-ICA datasets.\n', numel(files));

for i = 1:numel(files)
    try
        % --- Load pre-ICA dataset 
        EEG_for_ICA = pop_loadset('filename', files(i).name, 'filepath', inA);
        EEG_for_ICA = eeg_checkset(EEG_for_ICA, 'eventconsistency');
        fprintf('\n[ICA] %d/%d  %s\n', i, numel(files), EEG_for_ICA.setname);

        % --- Safety checks
        if any(~isfinite(EEG_for_ICA.data(:)))
            warning('Dataset %s contains NaN/Inf; skipping.', EEG_for_ICA.setname);
            continue
        end

        % --- Include all channels for ICA training
        chanind = 1:EEG_for_ICA.nbchan;

        % --- Effective rank (handles average-ref rank drop)
        if exist('eeg_rank','file')
            data_rank = eeg_rank(double(EEG_for_ICA.data(chanind,:)));
        else
            data_rank = rank(double(EEG_for_ICA.data(chanind,:)), 1e-7);
        end
        ncomps = min(numel(chanind), data_rank);
        if ncomps < 2
            warning('Rank < 2 for %s; skipping.', EEG_for_ICA.setname);
            continue
        end

        % --- Reproducible seed per subject 
        oldrng = rng; rng(100 + i, 'twister');

        % --- Run ICA: Infomax (extended), classic sphering (default), whole data
        tic;
        EEG_for_ICA = pop_runica(EEG_for_ICA, ...
            'icatype',  'runica', ...      % Infomax (gradient)
            'extended', 1, ...             % Extended ICA
            'chanind',  chanind, ...
            'pca',      ncomps, ...        % rank-aware
            'maxsteps', 512, ...           % # of ICA steps (BVA setting)
            'stop',     1e-7, ...          % convergence bound (BVA setting)
            'interrupt','off');            
        tmin = toc/60;
        rng(oldrng);

        % --- Bookkeeping
        if ~isfield(EEG_for_ICA,'etc') || ~isstruct(EEG_for_ICA.etc), EEG_for_ICA.etc = struct; end
        EEG_for_ICA.etc.ica = struct('chanind',chanind, 'rank',data_rank, ...
                                     'ncomps',size(EEG_for_ICA.icaweights,1), ...
                                     'maxsteps',512, 'stop',1e-7, 'extended',1, ...
                                     'sphering','on', 'runtime_min',tmin, ...
                                     'seed',100+i);
        fprintf('[ICA] %s: nIC=%d, rank=%d, time=%.1f min\n', ...
            EEG_for_ICA.setname, EEG_for_ICA.etc.ica.ncomps, data_rank, tmin);

        % --- Save dataset WITH ICA weights
        outName = [EEG_for_ICA.setname '_withICA.set'];
        pop_saveset(EEG_for_ICA, 'filename', outName, 'filepath', outB);

        % save info file
        info = EEG_for_ICA.etc.ica; %#ok<NASGU>
        save(fullfile(outB, [EEG_for_ICA.setname '_ica_info.mat']), 'info');

        % Free memory before next subject
        clear EEG_for_ICA

    catch ME
        warning('ICA failed for %s: %s', files(i).name, ME.message);
        continue
    end
end

fprintf('\n[Batch ICA] Done.\n');

%% 
 % % ---------------------------
 % % 10. Load ICA and identify ICA components
 % % ---------------------------
% Batch IC inspection: activations + ICLabel properties (manual review)

outB = '...\ICA_Results';
filesICA = dir(fullfile(outB, '*_withICA.set'));
fprintf('[Inspect] Found %d ICA datasets.\n', numel(filesICA));
if isempty(filesICA), error('No *_withICA.set found in %s', outB); end

set(0,'DefaultFigureWindowStyle','normal');

THR_art = 0.70;          % shortlist threshold (Eye/Muscle/Line/ChanNoise/Heart >= 0.70)

selDir = fullfile(outB, 'Manual_IC_Selections');
if ~isfolder(selDir), mkdir(selDir); end

for i = 1:numel(filesICA)
    fprintf('\n[Inspect %d/%d] %s\n', i, numel(filesICA), filesICA(i).name);


    % ---- 10.1. Load dataset with ICA
    EEG_ICA = pop_loadset('filename', filesICA(i).name, 'filepath', outB);
    EEG_ICA = eeg_checkset(EEG_ICA);
    % only run the following in batched processing, otherwise it will crash
    if isempty(EEG_ICA.icaweights)
        warning('No ICA weights in %s — skipping.', EEG_ICA.setname);
        continue
    end

    % number of components
    nIC = size(EEG_ICA.icaweights,1);

    % Ensure ICLabel exists
    if ~isfield(EEG_ICA,'etc') || ~isfield(EEG_ICA.etc,'ic_classification') ...
       || ~isfield(EEG_ICA.etc.ic_classification,'ICLabel')
        EEG_ICA = pop_iclabel(EEG_ICA, 'default');   % no auto-reject
    end


    % ---- 10.2. Build ICLabel-only shortlist (Eye/Muscle/Line/ChanNoise/Heart >= 0.70)
    
    classes = EEG_ICA.etc.ic_classification.ICLabel.classes;          
    P       = EEG_ICA.etc.ic_classification.ICLabel.classifications; 
    nIC     = size(EEG_ICA.icaweights,1);
    
    % robust class lookup
    normalize = @(s) regexprep(lower(s),'[\s_]+','');
    cls = cellfun(normalize, classes, 'uni', 0);
    findcls = @(name) find(strcmp(cls, normalize(name)), 1);
    
    iEye    = findcls('Eye');
    iMuscle = findcls('Muscle');
    iLine   = findcls('Line Noise');    if isempty(iLine),   iLine   = findcls('LineNoise');    end
    iChanN  = findcls('Channel Noise'); if isempty(iChanN),  iChanN  = findcls('ChannelNoise'); end
    iHeart  = findcls('Heart');         if isempty(iHeart),  iHeart  = findcls('Heartbeat');    end
    
    colEye   = P(:,iEye);
    colMus   = P(:,iMuscle);
    colLine  = P(:,iLine);
    colChanN = P(:,iChanN);
    colHeart = zeros(nIC,1); if ~isempty(iHeart), colHeart = P(:,iHeart); end
    
    artifactScore = max([colEye colMus colLine colChanN colHeart], [], 2, 'omitnan');
    
    % shortlist = ICs meeting threshold, **sorted by IC number (ascending)**
    shortlist = find(artifactScore >= THR_art);
    shortlist = sort(shortlist);
    

   % ---- 10.03 Show shortlist in its own window (IC as integers, no row numbers)
    try, close(hShort); end %#ok<TRYNC>
    
    if ~isempty(shortlist)
        icStr  = arrayfun(@(k) sprintf('%d',k), shortlist, 'uni', false);  % IC as integer strings
        tblData = [ icStr, ...
                    num2cell(round(colEye(shortlist),  2)), ...
                    num2cell(round(colMus(shortlist),  2)), ...
                    num2cell(round(colLine(shortlist), 2)), ...
                    num2cell(round(colChanN(shortlist),2)), ...
                    num2cell(round(colHeart(shortlist),2)), ...
                    num2cell(round(artifactScore(shortlist),2)) ];
    
        hShort = figure('Name', sprintf('%s | Shortlist (>= %.0f%%)', EEG_ICA.setname, THR_art*100), ...
                        'NumberTitle','off', 'MenuBar','none', 'ToolBar','none', ...
                        'Position', [100 100 600 320]);
        uitable('Parent', hShort, ...
                'Data', tblData, ...
                'ColumnName', {'IC','Eye','Muscle','Line','ChanNoise','Heart','Max'}, ...
                'ColumnEditable', false(1,7), ...
                'RowName', [], ...                % <-- hide row numbers
                'Units','normalized', 'Position',[0 0 1 1]);
        try, movegui(hShort, 'east'); end
    else
        hShort = figure('Name', sprintf('%s | Shortlist (none >= %.0f%%)', EEG_ICA.setname, THR_art*100), ...
                        'NumberTitle','off', 'MenuBar','none', 'ToolBar','none', ...
                        'Position', [100 100 420 80]);
        uicontrol('Style','text','Parent',hShort,'Units','normalized','Position',[.05 .2 .9 .6], ...
                  'String','No components met the shortlist threshold.', 'FontSize', 11);
        try, movegui(hShort, 'east'); end
    end
    

    % ---- 10.04 IC activations viewer 
    ica_act  = eeg_getica(EEG_ICA);
    iclabels = arrayfun(@(k) sprintf('IC%03d',k), 1:nIC, 'UniformOutput', false);
    iclocs   = struct('labels', iclabels);
    actTitle = sprintf('%s | IC activations (12/page)', EEG_ICA.setname);
    eegplot(ica_act, 'srate',EEG_ICA.srate, 'winlength',8, 'dispchans',12, 'spacing',15, ...
            'eloc_file',iclocs, 'events',EEG_ICA.event, 'title', actTitle);
    drawnow; pause(0.05);
    hAct = findall(0,'type','figure','tag','EEGPLOT','Name',actTitle);
    if ~isempty(hAct)
        hChk = findall(hAct,'style','checkbox');
        for k = 1:numel(hChk)
            txt = get(hChk(k),'String');
            if ischar(txt) && contains(lower(txt),'remove dc','IgnoreCase',true)
                set(hChk(k),'Value',1);
                cb = get(hChk(k),'Callback');
                if isa(cb,'function_handle'), cb(hChk(k),[]);
                elseif ischar(cb), eval(cb);
                end
                break
            end
        end
        try, movegui(hAct,'west'); end
    end
    
    % ---- 10.5 Pre-mark shortlist BEFORE registering the dataset
    nIC = size(EEG_ICA.icaweights,1);
    if ~isfield(EEG_ICA,'reject') || ~isstruct(EEG_ICA.reject) || ~isfield(EEG_ICA.reject,'gcompreject')
        EEG_ICA.reject.gcompreject = false(1,nIC);
    else
        EEG_ICA.reject.gcompreject(:) = false;
    end
    if exist('shortlist','var') && ~isempty(shortlist)
        EEG_ICA.reject.gcompreject(shortlist) = true;
    end
    
    % ---- Register dataset in EEGLAB, sync base, open selector once
    eeglab;
    [ALLEEG, EEG_ICA, CURRENTSET] = eeg_store([], EEG_ICA, 0);
    assignin('base','EEG', EEG_ICA);                 
    
    % One-window selector (built-in paging via Next/Prev)
    [EEG_ICA, ~] = pop_selectcomps(EEG_ICA, 1:nIC);  
    hSel = findall(0,'type','figure','-regexp','Name','Component.*selection|IC.*selection');
    if isempty(hSel), hSel = gcf; end
    if isgraphics(hSel(1)), waitfor(hSel(1)); end
    
    % --- AFTER the selector window is closed, pull the edited dataset back
    EEG_from_gui = [];
    try
        EEG_from_gui = evalin('base','EEG');  
    catch
    end
    if isstruct(EEG_from_gui) && isfield(EEG_from_gui,'reject') && ...
            isfield(EEG_from_gui.reject,'gcompreject') && ...
            numel(EEG_from_gui.reject.gcompreject) == nIC
        EEG_ICA = EEG_from_gui;                % take the GUI-edited struct
    else
        % fallback: EEGLAB container copy
        EEG_ICA = ALLEEG(CURRENTSET);
    end
    
    % sanity print
    fprintf('Marked for reject (now): %s\n', mat2str(find(EEG_ICA.reject.gcompreject)));
    
    % --- Now persist back into EEGLAB container (overwrite the set)
    ALLEEG(CURRENTSET) = EEG_ICA;
    assignin('base','EEG', EEG_ICA);   % keep base in sync for safety

    
    % ---- 10.6. Report rejected ICs
    rejIC = find(EEG_ICA.reject.gcompreject);
    fprintf('[%s] Rejected ICs: %s\n', EEG_ICA.setname, mat2str(rejIC));
    
    % ---- 10.7. Save a copy with the marks + verify
    baseSafe = regexprep(EEG_ICA.setname,'[<>:"/\\|?*]+','_');
    pop_saveset(EEG_ICA, 'filename', [baseSafe '_withICA_marks.set'], 'filepath', selDir);
    EEG_chk = pop_loadset('filename', [baseSafe '_withICA_marks.set'], 'filepath', selDir);
    fprintf('[%s] Rejected ICs (saved file): %s\n', EEG_chk.setname, mat2str(find(EEG_chk.reject.gcompreject)));
    
    % Tidy windows 
    if exist('hAct','var')
        try
            hValid = hAct(isgraphics(hAct));
            if ~isempty(hValid), close(hValid); end
        end
    end
    if exist('hShort','var')
        try
            if isgraphics(hShort), close(hShort); end
        end
    end
    if exist('hPreviewAll','var')
        try
            if isgraphics(hPreviewAll), close(hPreviewAll); end
        end
    end
    if exist('hPreviewLater','var')
        try
            if isgraphics(hPreviewLater), close(hPreviewLater); end
        end
    end
end

fprintf('\n[Inspect] Done. Selections saved to: %s\n', selDir);

%% 
 % % ---------------------------
 % % 11. Remove ICA components from EEG (0.1 Hz filtered) data
 % % ---------------------------

inAnalysis = '...\Pre_ICA_Sets';
inMarks    = '...\Manual_IC_Selections';
outPostICA = 'S:\KJP_Neurophys\Studien\STIPED\WP5 Querschnitt (SoKoASS)\6_Projekte\1_Posner\2_Paper Leonie - oscillations\Auswertung\Post_ICA_Sets';
if ~isfolder(outPostICA), mkdir(outPostICA); end

markFiles = dir(fullfile(inMarks, '*_withICA_marks.set'));
fprintf('[Apply ICA] Found %d marked ICA sets.\n', numel(markFiles));

for i = 1:numel(markFiles)
    try
        % --- 1) Load marked ICA dataset (has weights + gcompreject)
        EEG_ICA = pop_loadset('filename', markFiles(i).name, 'filepath', inMarks);
        EEG_ICA = eeg_checkset(EEG_ICA, 'eventconsistency');
        if isempty(EEG_ICA.icaweights)
            warning('No ICA weights in %s; skipping.', EEG_ICA.setname);
            continue
        end
        rejIC = find(EEG_ICA.reject.gcompreject);
        if isempty(rejIC)
            fprintf('[%s] No ICs marked; will still transfer ICA.\n', EEG_ICA.setname);
        else
            fprintf('[%s] Removing ICs: %s\n', EEG_ICA.setname, mat2str(rejIC));
        end

        % --- 2) Locate matching 0.1 Hz analysis set
        subj = EEG_ICA.subject;
        if isempty(subj)
            % fallback: derive from setname (strip trailing "_withICA_marks")
            subj = regexprep(EEG_ICA.setname, '_withICA.*$', '');
        end
        cand = fullfile(inAnalysis, [subj '_preICA_analysis.set']);
        if ~isfile(cand)
            % fuzzy search: any analysis file that contains the subject string
            allA = dir(fullfile(inAnalysis, '*_preICA_analysis.set'));
            hits = allA(contains({allA.name}, subj, 'IgnoreCase', true));
            if isempty(hits)
                warning('[%s] Could not find matching analysis set for subject "%s".', EEG_ICA.setname, subj);
                continue
            elseif numel(hits) > 1
                % pick the closest by longest common substring with setname
                [~, pick] = max(cellfun(@(n) numel(intersect(lower(n), lower(EEG_ICA.setname))), {hits.name}));
                cand = fullfile(inAnalysis, hits(pick).name);
            else
                cand = fullfile(inAnalysis, hits(1).name);
            end
        end
        fprintf('  ↳ Analysis set: %s\n', cand);

        % Robust load with filepath split (Windows-safe)
        [candPath, candFile, candExt] = fileparts(cand);
        EEG = pop_loadset('filename', [candFile candExt], 'filepath', candPath);
        EEG = eeg_checkset(EEG);

        % --- 3) Make sure channel sets match; reorder analysis channels if needed
        labsA   = string({EEG.chanlocs.labels});
        labsICA = string({EEG_ICA.chanlocs.labels});
        if ~isequal(labsA, labsICA)
            % same set, different order?
            [tf, idx] = ismember(labsICA, labsA);
            if all(tf)
                EEG = pop_select(EEG, 'channel', cellstr(labsICA));  % reorder to match ICA
                labsA = string({EEG.chanlocs.labels});
            else
                warning('[%s] Channel sets differ; cannot transfer ICA safely. Skipping.', EEG_ICA.setname);
                continue
            end
        end

        % --- 4) Transfer ICA decomposition to analysis EEG
        EEG.icawinv     = EEG_ICA.icawinv;
        EEG.icaweights  = EEG_ICA.icaweights;
        EEG.icasphere   = EEG_ICA.icasphere;
        if isfield(EEG_ICA, 'icachansind'), EEG.icachansind = EEG_ICA.icachansind; end
        EEG = eeg_checkset(EEG);

        % --- 5) Remove marked components (if any)
        if ~isempty(rejIC)
            EEG_clean = pop_subcomp(EEG, rejIC, 0);   % 0 = remove ICs from data
        else
            EEG_clean = EEG;
        end
        EEG_clean = eeg_checkset(EEG_clean);

        % keep a small audit trail
        if ~isfield(EEG_clean,'etc') || ~isstruct(EEG_clean.etc), EEG_clean.etc = struct; end
        EEG_clean.etc.ic_apply = struct( ...
            'source_set', EEG_ICA.setname, ...
            'subject',    subj, ...
            'removedICs', rejIC, ...
            'nRemoved',   numel(rejIC), ...
            'timestamp',  datestr(now, 'yyyy-mm-dd HH:MM:SS'));

        % Sanity check: Which ICs were removed?
        if isfield(EEG_clean,'etc') && isfield(EEG_clean.etc,'ic_apply')
            disp(EEG_clean.etc.ic_apply)
        end
        % ICA matrices present?
        whos EEG.icawinv EEG.icaweights EEG.icasphere

        % --- 6) Save post-ICA analysis set
        outName = sprintf('%s_postICA_analysis.set', regexprep(subj,'[<>:"/\\|?*]+','_'));
        pop_saveset(EEG_clean, 'filename', outName, 'filepath', outPostICA);
        fprintf(' Saved: %s\n', fullfile(outPostICA, outName));


        % free memory for next subject
        clear EEG EEG_clean EEG_ICA

    catch ME
        warning('Apply failed for %s: %s', markFiles(i).name, ME.message);
        continue
    end
end

fprintf('\n[Apply ICA] Done. Output in: %s\n', outPostICA);


%% 
 % % ---------------------------
 % % 12.  Interpolation + Epoching + Baseline + Trial rejection
 % % ---------------------------

inPostICA  = '...\Post_ICA_Sets';
outInterp  = '...\Post_ICA_Interp_FCz';
outEpochs  = '...\Epochs';
if ~isfolder(outInterp), mkdir(outInterp); end
if ~isfolder(outEpochs), mkdir(outEpochs); end

files = dir(fullfile(inPostICA, '*_postICA_analysis.set'));
fprintf('[Stage D/E] Found %d post-ICA analysis sets.\n', numel(files));
set(0,'DefaultFigureWindowStyle','normal');

for i = 1:numel(files)
    try
        % ---------- Load post-ICA set ----------
        EEG = pop_loadset('filename', files(i).name, 'filepath', inPostICA);
        EEG = eeg_checkset(EEG);

        subj = EEG.subject; if isempty(subj), subj = EEG.setname; end
        safeSubj = regexprep(subj,'[<>:"/\\|?*]+','_');

        % --- Common: load full montage once -> as channel locations for
        % rejected channels have been deleted
        chanlocs_file = '...\64_extended20-10.elp';
        fullLocs  = readlocs(chanlocs_file);
        labsFull  = {fullLocs.labels};

        isIO_full  = strcmpi(labsFull,'IO');
        isFCz_full = strcmpi(labsFull,'FCz');

        tmplScalpLocs = fullLocs(~(isIO_full | isFCz_full));      % scalp-only template
        tmplScalpLabs = {tmplScalpLocs.labels};
        
        % -------- Current dataset labels --------
        presentLabs   = {EEG.chanlocs.labels};
        presentScalp  = presentLabs(~strcmpi(presentLabs,'IO') & ~strcmpi(presentLabs,'FCz'));

        % --- Drop ICA fields 
        EEG.icaweights  = []; EEG.icasphere = []; EEG.icawinv = []; EEG.icaact = [];
        if isfield(EEG,'icachansind'), EEG.icachansind = []; end
        
        nb_before = EEG.nbchan;

        % --- BEFORE interpolation: capture current labels
        labs_before = {EEG.chanlocs.labels};
        
        % -------- Find missing scalp channels (ignore IO & FCz) --------
        % AFTER (case-insensitive):
        tmplScalpLabs_s = string(tmplScalpLabs);
        presentScalp_s  = string(presentScalp);
        missingMask = ~ismember(lower(tmplScalpLabs_s), lower(presentScalp_s));
        missingLabs = cellstr(tmplScalpLabs_s(missingMask));
        
        if ~isempty(missingLabs)
            % only interpolate the actually missing scalp channels
            missLocs = tmplScalpLocs(ismember(tmplScalpLabs, missingLabs));
            fprintf('  Interpolating missing scalp channels (%d): %s\n', numel(missingLabs), strjoin(missingLabs, ', '));
            EEG = pop_interp(EEG, missLocs, 'spherical');          % pass a *subset* chanlocs struct → EEGLAB adds just these
            EEG = eeg_checkset(EEG);
        else
            fprintf('  No scalp channels missing — skipping interpolation.\n');
        end
        
        % documentation
        EEG.etc.interpMissing = missingLabs;
        EEG.etc.interpMethod  = 'spherical';

        
        % -------- Average re-ref (exclude IO), and add FCz via refloc  --------
        iFCz_full = find(strcmpi(labsFull,'FCz'), 1);
        fcZloc    = fullLocs(iFCz_full);
        if ~isfield(fcZloc,'ref'),    fcZloc.ref = [];    end
        if ~isfield(fcZloc,'urchan'), fcZloc.urchan = []; end

        % Before reref: compute mean across the channels that will be included in the average
        use = find(~ismember({EEG.chanlocs.labels}, {'IO','FCz'}));   
        refSig_pre = mean(EEG.data(use,:), 1);   
        
        % keep IO if present, but exclude it from the average
        ixIO = find(strcmpi({EEG.chanlocs.labels}, 'IO'), 1);
        if isempty(ixIO)
            EEG = pop_reref(EEG, [], 'refloc', fcZloc);
        else
            EEG = pop_reref(EEG, [], 'exclude', ixIO, 'refloc', fcZloc);
        end
        EEG = eeg_checkset(EEG);

        % enforce a single final channel order across all files
        % template scalp order (no IO, no FCz)
        targetScalpOrder = tmplScalpLabs;  % from your template file
        
        % build final order: scalp (template order) -> FCz -> (optional IO last)
        hasIO = ~isempty(find(strcmpi({EEG.chanlocs.labels}, 'IO'), 1));
        finalOrder = [targetScalpOrder, {'FCz'}];
        if hasIO
            finalOrder = [finalOrder, {'IO'}];
        end
        
        % apply exact order (and drop any extras not in template)
        EEG = pop_select(EEG, 'channel', finalOrder);
        EEG = eeg_checkset(EEG);

                    
        % --- POST checks
        fprintf('Re-ref done. Channels before: %d, after: %d\n', nb_before, EEG.nbchan);

        % After reref: FCz should equal V_FCz - V_avg = -(V_avg - V_FCz) = -refSig_pre
        iFCz = find(strcmpi({EEG.chanlocs.labels},'FCz'),1);
        rmsDiff_plus  = rms(double(EEG.data(iFCz,:) -  refSig_pre));
        rmsDiff_minus = rms(double(EEG.data(iFCz,:) +  refSig_pre));
        fprintf('RMS(FCz - mean_pre)  = %.4g µV\n', rmsDiff_plus);
        fprintf('RMS(FCz + mean_pre)  = %.4g µV (this one should be ~0)\n', rmsDiff_minus);

        % FCz presence & near-zero data
        iFCz_now = find(strcmpi({EEG.chanlocs.labels},'FCz'),1);
        if isempty(iFCz_now)
            warning('FCz missing after re-ref!');
        else
            maxFCz = max(abs(EEG.data(iFCz_now,:)));
            fprintf('  FCz index %d | max|signal| = %.3g µV \n', iFCz_now, maxFCz);
        end

        % Final montage consistency: after interpolate+retrieve, data should have
        % exactly the full montage (template) channels present 
        missing_from_template = setdiff(labsFull, {EEG.chanlocs.labels});
        extra_not_in_template = setdiff({EEG.chanlocs.labels}, labsFull);
        fprintf('  Missing from template: %s\n', strjoin(missing_from_template, ', '));
        fprintf('  Extra not in template: %s\n', strjoin(extra_not_in_template, ', '));
        
        % (Optional) assert expected final count equals template count
        if isempty(missing_from_template) && isempty(extra_not_in_template)
            assert(EEG.nbchan == numel(fullLocs), 'Channel count mismatch after re-ref.');
        end

        % save
        outNameCont = sprintf('%s_postICA_interpFCz.set', safeSubj);
        pop_saveset(EEG, 'filename', outNameCont, 'filepath', outInterp);
        fprintf('  Saved: %s\n', fullfile(outInterp, outNameCont));
    

        %% ===================== SEGMENTATION (TF-focused) =====================
        % Includes:
        %  - define event codes & groups
        %  - build a trial table on the continuous data (cue->target pairing)
        %  - add behavior (resp side, RT, accuracy, error type) to that table
        %  - decide epoch spans with enough padding for TF
        %  - epoch cue-locked and target-locked data
        %  - attach per-epoch metadata including TF baseline windows (same physiology baseline)
        
        % --- Event codes  ---
        CUE_L_short = {'110'};   CUE_L_long = {'210'};
        CUE_R_short = {'120'};   CUE_R_long = {'220'};
        CUE_L_CATCH = {'10'};    CUE_R_CATCH = {'20'};
        
        TARG_L_valid_long    = {'211'};  TARG_L_valid_short   = {'111'};
        TARG_L_invalid_long  = {'221'};  TARG_L_invalid_short = {'121'};
        TARG_R_valid_long    = {'222'};  TARG_R_valid_short   = {'122'};
        TARG_R_invalid_long  = {'212'};  TARG_R_invalid_short = {'112'};
        
        % Uncued targets 
        TARG_L_uncued = {'1'};   TARG_R_uncued = {'2'};
        
        % --- Groups (convenience sets) ---
        CUE_L_ALL = [CUE_L_short CUE_L_long CUE_L_CATCH];
        CUE_R_ALL = [CUE_R_short CUE_R_long CUE_R_CATCH];
        CUE_ALL   = [CUE_L_ALL CUE_R_ALL];
        
        TARG_L_ALL = [TARG_L_valid_long TARG_L_valid_short TARG_L_invalid_long TARG_L_invalid_short TARG_L_uncued];
        TARG_R_ALL = [TARG_R_valid_long TARG_R_valid_short TARG_R_invalid_long TARG_R_invalid_short TARG_R_uncued];
        TARG_ALL   = [TARG_L_ALL TARG_R_ALL];
        
        % --- Responses (for behavior metrics that will be attached to the trial table) ---
        RESP_L = {'7'};        % left button code
        RESP_R = {'8'};        % right button code
        rt_min_ms = 200;       % anticipations faster than this
        rt_max_ms = 1000;      % too-late responses
        fa_window_ms = 1000;   % catch trials: look this long after the cue for a false alarm
        
        % --- Pairing window (cue → next target) on the continuous stream ---
        pair_window_sec = 0.5;  % covers variable SOA with margin
        
        % --- TF padding & guard ---
        % Padding protects TF from edge artifacts at lowest frequency
        % (alpha)
        % Guard keeps target activity from smearing into cue-locked analysis.
        f_min = 8; 
        n_cycles_low = 4;                   
        pad_ms = ceil(1000 * n_cycles_low / f_min);   % ~500 ms padding
        guard_ms = 200;                               % keep-out before target (cue-locked)
        
        % --- Epoch spans (with padding) ---
        % Cue: earlier [-1.0, +1.5] already includes baseline & generous right pad.
        cue_epoch_sec  = [-(0.5 + pad_ms/1000),  1.0 + pad_ms/1000];  % -> ~[-1.0, +1.5] s
        
        % Target: start far enough left to include the mapped pre-cue baseline
        % [-SOA-500, -SOA-200] PLUS pad. Compute SOA_max from the trial table below
        % and then set targ_epoch_sec dynamically
        
        % --- TF baseline window  ---
        % Baseline will be used in TF (dB) later
        % -> segment baseline (including padding) to target (including padding), 
        % compute TF over the whole segment, and then only use the time windows 
        % of interest for dB computation
        tf_base_precue_ms = [-500 -200];

        % --- Normalize event types to numeric strings (e.g., 'S222' -> '222', 'S 8' -> '8')
        % print a normalized inventory without changing EEG.event.type
        keys_raw  = string({EEG.event.type});
        keys_norm = regexprep(upper(keys_raw), '\D', '');   % strip non-digits: 'S222'->'222', 'S 8'->'8'
        emptyMask = (keys_norm == "");
        keys_norm(emptyMask) = keys_raw(emptyMask);         % keep non-numeric labels (e.g., 'BOUNDARY')

        [u,~,ic] = unique(keys_norm);
        counts   = accumarray(ic,1);
        fprintf('[Events] %d unique types after normalization-like view:\n', numel(u));
        for k = 1:numel(u)
            fprintf('  %4d  %s\n', counts(k), u{k});
        end
        
        EEG = eeg_checkset(EEG, 'eventconsistency');

    
        
        % ---------------- Build / reuse trial table on the continuous EEG (incl. UNCUED trials) ----------------
        % Helper: normalize event labels by keeping only digits (so 'S 7' -> '7', 'S222' -> '222')
        toKey = @(x) regexprep(upper(string(x)),'\D','');
        on    = @(t,arr) ismember(toKey(t), toKey(arr));
        
        % Define cued-target set for pairing, and uncued explicitly
        TARG_CUED = [TARG_L_valid_long TARG_L_valid_short ...
                     TARG_L_invalid_long TARG_L_invalid_short ...
                     TARG_R_valid_long TARG_R_valid_short ...
                     TARG_R_invalid_long TARG_R_invalid_short];
        TARG_UNCUED = [TARG_L_uncued TARG_R_uncued];
        
        % (Re)build the trial table (force rebuild to ensure uncued targets are included)
        ev = EEG.event; nE = numel(ev); sr = EEG.srate;
        pair_window_pts = round(pair_window_sec * sr);
        
        TT = struct('cue_ev',{},'cue_lat',{},'cue_side',{},'cue_soa',{}, ...
                    'cue_urev',{}, ...
                    'targ_ev',{},'targ_lat',{},'targ_side',{},'targ_soa',{},'validity',{}, ...
                    'targ_urev',{});
        
        usedTargs = false(1, nE);  % mark targets already paired to a cue

        
        % ---------- Pass 1: one trial per CUE, pairing only to a cued target within the window ----------
        e = 1;
        while e <= nE
            typ = ev(e).type;
            if on(typ, CUE_ALL)
                % Start a new cue-based trial
                rec.cue_ev  = e;
                rec.cue_lat = round(ev(e).latency);
        
                % cue_side: L/R/N
                if on(typ, [CUE_L_short CUE_L_long CUE_L_CATCH])
                    rec.cue_side = "L";
                elseif on(typ, [CUE_R_short CUE_R_long CUE_R_CATCH])
                    rec.cue_side = "R";
                else
                    rec.cue_side = "N";
                end
        
                % cue_soa: short/long/catch/""
                if on(typ, [CUE_L_short CUE_R_short])
                    rec.cue_soa = "short";
                elseif on(typ, [CUE_L_long CUE_R_long])
                    rec.cue_soa = "long";
                elseif on(typ, [CUE_L_CATCH CUE_R_CATCH])
                    rec.cue_soa = "catch";
                else
                    rec.cue_soa = "";
                end
        
                % cue_urev (safe)
                if isfield(ev,'urevent') && ~isempty(ev(e).urevent)
                    rec.cue_urev = ev(e).urevent;
                else
                    rec.cue_urev = NaN;
                end
        
                % Init target-related fields
                rec.targ_ev   = [];
                rec.targ_lat  = NaN;
                rec.targ_side = "N";
                rec.targ_soa  = "";
                rec.validity  = "cue_only";
                rec.targ_urev = NaN;
        
                % Find the first CUED target within the pairing window
                j = e + 1;
                while j <= nE && round(ev(j).latency) <= rec.cue_lat + pair_window_pts
                    if on(ev(j).type, TARG_CUED)
                        rec.targ_ev  = j;
                        rec.targ_lat = round(ev(j).latency);
        
                        % targ_side: L / R
                        if on(ev(j).type, [TARG_L_valid_long TARG_L_valid_short TARG_L_invalid_long TARG_L_invalid_short])
                            rec.targ_side = "L";
                        else
                            rec.targ_side = "R";
                        end
        
                        % targ_soa: short / long
                        if on(ev(j).type, [TARG_L_valid_short TARG_L_invalid_short TARG_R_valid_short TARG_R_invalid_short])
                            rec.targ_soa = "short";
                        else
                            rec.targ_soa = "long";
                        end
        
                        % validity: valid / invalid
                        if on(ev(j).type, [TARG_L_valid_long TARG_L_valid_short TARG_R_valid_long TARG_R_valid_short])
                            rec.validity = "valid";
                        else
                            rec.validity = "invalid";
                        end
        
                        % targ_urev (safe)
                        if isfield(ev,'urevent') && ~isempty(ev(j).urevent)
                            rec.targ_urev = ev(j).urevent;
                        else
                            rec.targ_urev = NaN;
                        end
        
                        usedTargs(j) = true;  % mark this target as used
                        break
                    end
                    j = j + 1;
                end
        
                TT(end+1) = rec; %#ok<AGROW>
            end
            e = e + 1;
        end
        
        % ---------- Pass 2: add one trial per *UNCUED* target (codes 1/2) ----------
        nAddedUncued = 0;
        for j = 1:nE
            if usedTargs(j), continue; end
            if on(ev(j).type, TARG_UNCUED)
                rec.cue_ev  = [];                 % no cue
                rec.cue_lat = NaN;
                rec.cue_side= "N";
                rec.cue_soa = "";
                rec.cue_urev= NaN;
        
                rec.targ_ev  = j;
                rec.targ_lat = round(ev(j).latency);
        
                if on(ev(j).type, TARG_L_uncued), rec.targ_side = "L"; else, rec.targ_side = "R"; end
                rec.targ_soa = "";                 % uncued → no short/long
                rec.validity = "no_cue";
        
                if isfield(ev,'urevent') && ~isempty(ev(j).urevent)
                    rec.targ_urev = ev(j).urevent;
                else
                    rec.targ_urev = NaN;
                end
        
                TT(end+1) = rec; %#ok<AGROW>
                nAddedUncued = nAddedUncued + 1;
            end
        end
        
        % ---------- Sort all trials by their onset time in the continuous stream ----------
        % (cue time if present, otherwise target time)
        onset_pts = zeros(1, numel(TT));
        for k = 1:numel(TT)
            if ~isnan(TT(k).cue_lat) && ~isempty(TT(k).cue_lat)
                onset_pts(k) = TT(k).cue_lat;
            else
                onset_pts(k) = TT(k).targ_lat;
            end
        end
        [~, ord] = sort(onset_pts, 'ascend');
        TT = TT(ord);
        
        % Write back
        if ~isfield(EEG,'etc') || ~isstruct(EEG.etc), EEG.etc = struct; end
        EEG.etc.trialtable = TT;
        
        fprintf('[Trials] Built trial table: %d cue-based + %d uncued-target trials (total=%d).\n', ...
                sum(~cellfun(@isempty,{TT.cue_ev})), nAddedUncued, numel(TT));

        
        
        % ---------------- Add behavior to the trial table (RT, accuracy, errors) --------
        TT = EEG.etc.trialtable;
        sr = EEG.srate;
        rt_min_pts = round(rt_min_ms/1000 * sr);
        rt_max_pts = round(rt_max_ms/1000 * sr);
        fa_win_pts = round(fa_window_ms/1000 * sr);
        
        for ii = 1:numel(TT)
            % init fields
            TT(ii).resp_ev    = [];
            TT(ii).resp_lat   = NaN;
            TT(ii).resp_side  = "";
            TT(ii).resp_urev  = NaN;
            TT(ii).rt_ms      = NaN;
            TT(ii).acc        = false;
            TT(ii).error_type = "";
        
            if ~isempty(TT(ii).targ_ev)  % ---------- target-present trial ----------
                tlat = TT(ii).targ_lat;
                j = TT(ii).targ_ev + 1; first = [];
                while j <= numel(EEG.event) && round(EEG.event(j).latency) <= (tlat + rt_max_pts)
                    if on(EEG.event(j).type, [RESP_L RESP_R]), first = j; break; end
                    j = j + 1;
                end
                if isempty(first)
                    TT(ii).error_type = "miss";  % no response in window
                else
                    TT(ii).resp_ev   = first;
                    TT(ii).resp_lat  = round(EEG.event(first).latency);
                    if isfield(EEG.event,'urevent') && ~isempty(EEG.event(first).urevent)
                        TT(ii).resp_urev = EEG.event(first).urevent;
                    end
                    % side mapping (no string arithmetic)
                    if on(EEG.event(first).type, RESP_L)
                        TT(ii).resp_side = "L";
                    elseif on(EEG.event(first).type, RESP_R)
                        TT(ii).resp_side = "R";
                    else
                        TT(ii).resp_side = "";
                    end
                    % RT + classification
                    rt_pts       = TT(ii).resp_lat - tlat;
                    TT(ii).rt_ms = 1000 * rt_pts / sr;
                    if rt_pts < rt_min_pts
                        TT(ii).error_type = "anticipation";
                    elseif rt_pts > rt_max_pts
                        TT(ii).error_type = "late";
                    else
                        TT(ii).acc = (TT(ii).resp_side == TT(ii).targ_side);
                        if TT(ii).acc
                            TT(ii).error_type = "correct";
                        else
                            TT(ii).error_type = "wrong_side";
                        end
                    end
                end
        
            else                  % ---------- catch trial (no target) ----------
                j = TT(ii).cue_ev + 1; fa = false; first = [];
                while j <= numel(EEG.event) && round(EEG.event(j).latency) <= (TT(ii).cue_lat + fa_win_pts)
                    if on(EEG.event(j).type, [RESP_L RESP_R]), fa = true; first = j; break; end
                    j = j + 1;
                end
                if fa
                    TT(ii).resp_ev   = first;
                    TT(ii).resp_lat  = round(EEG.event(first).latency);
                    if isfield(EEG.event,'urevent') && ~isempty(EEG.event(first).urevent)
                        TT(ii).resp_urev = EEG.event(first).urevent;
                    end
                    if on(EEG.event(first).type, RESP_L)
                        TT(ii).resp_side = "L";
                    elseif on(EEG.event(first).type, RESP_R)
                        TT(ii).resp_side = "R";
                    else
                        TT(ii).resp_side = "";
                    end
                    TT(ii).acc = false; TT(ii).error_type = "false_alarm";
                else
                    TT(ii).acc = true;  TT(ii).error_type = "correct_reject";
                end
            end
        end
        EEG.etc.trialtable = TT;    % write back

        % -------- Per-trial CSV (audit) --------
        try
            nTT = numel(TT);
            T = struct2table(TT);
        
            % Add sequential trial index in continuous order
            T.trial = (1:nTT).';
            T = movevars(T, 'trial', 'Before', 1);
        
            % Make sure string-ish fields are strings (avoids cellstr/char mix in CSV)
            for vn = string(T.Properties.VariableNames)
                if iscell(T.(vn)) && ~isempty(T.(vn)) && ischar(T.(vn){1})
                    T.(vn) = string(T.(vn));
                end
            end
        
            if ~exist('outEpochs','var') || isempty(outEpochs), outEpochs = pwd; end
            subj = EEG.subject; if isempty(subj), subj = EEG.setname; end
            safeSubj = regexprep(subj,'[<>:"/\\|?*]+','_');
            outCSV = fullfile(outEpochs, [safeSubj '_trialtable.csv']);
        
            writetable(T, outCSV, 'Delimiter',';');  % German Excel friendly
            fprintf('  Wrote per-trial CSV: %s (n=%d trials)\n', outCSV, nTT);
        catch ME
            warning('Could not write per-trial CSV: %s', ME.message);
        end
        
        
        % ---------------- Decide target epoch span dynamically (covers baseline+pad) ----
        hasTarg = ~arrayfun(@(x) isempty(x.targ_ev), TT);
        if any(hasTarg)
            d_pts = [TT(hasTarg).targ_lat] - [TT(hasTarg).cue_lat];
            SOA_ms_all = 1000 * (d_pts / sr);
            SOA_max_ms = max(SOA_ms_all);                 
        else
            SOA_max_ms = 200;                              % fallback if no targets found
        end
        
        % Left side must reach: -(SOA_max + 500) ms (baseline start) minus pad
        targ_left_sec  = -(SOA_max_ms + 500 + pad_ms)/1000;   % e.g., -(200+500+500)/1000 = -1.2 s
        % Right side must cover analysis (+800 ms) plus pad
        targ_right_sec = (800 + pad_ms)/1000;                 % 0.8 + 0.5 = 1.3 s
        targ_epoch_sec = [targ_left_sec, targ_right_sec];
        
       %% --- Preflight before epoching  ---

        % Normalize event type once to bare numeric (digits) and keep an original backup
        if ~isfield(EEG,'etc') || ~isfield(EEG.etc,'event_type_original')
            EEG.etc.event_type_original = {EEG.event.type};
            for ii = 1:numel(EEG.event)
                s = upper(strtrim(char(string(EEG.event(ii).type))));
                tok = regexp(s,'\d+','match');    % grab digit blocks
                if ~isempty(tok)
                    EEG.event(ii).type = tok{end}; % keep last numeric block, e.g., '222'
                else
                    EEG.event(ii).type = s;        % non-numeric (e.g., 'BOUNDARY')
                end
            end
            EEG = eeg_checkset(EEG,'eventconsistency');
        end
        
        % Helpers for normalized membership tests
        toKey = @(x) regexprep(upper(string(x)),'\D','');  % keep only digits
        keys  = toKey({EEG.event.type});                   % normalized view of the stream
        
        % 1) How many targets are present after normalization?
        isTarg  = ismember(keys, toKey(TARG_ALL));
        targIdx = find(isTarg);
        fprintf('[Preflight] Found %d target events matching TARG_ALL.\n', numel(targIdx));
        if isempty(targIdx)
            error('No target events match TARG_ALL. Check your code lists vs. EEG.event.type.');
        end
        
        % 2) Which targets are feasible given the target epoch window? (EXACT test)
        sr           = EEG.srate;
        targ_left_sec  = targ_epoch_sec(1);   % negative
        targ_right_sec = targ_epoch_sec(2);   % positive
        left_samp_t    = targ_left_sec  * sr; % negative
        right_samp_t   = targ_right_sec * sr; % positive
        
        evLatExactT = double([EEG.event(targIdx).latency]);   % 1-based, fractional allowed
        okT_exact   = (evLatExactT + left_samp_t  >= 1) & (evLatExactT + right_samp_t <= EEG.pnts);
        
        % (optional) rounded version just for logging
        tmin_samp_r = 1 - round(targ_left_sec  * sr);
        tmax_samp_r = EEG.pnts - round(targ_right_sec * sr);
        evLatRoundT = round(evLatExactT);
        okT_round   = (evLatRoundT >= tmin_samp_r) & (evLatRoundT <= tmax_samp_r);
        
        fprintf('[Preflight] Targets feasible (exact) with window [%.3f %.3f] s — feasible: %d / %d (rounded=%d)\n', ...
                targ_left_sec, targ_right_sec, sum(okT_exact), numel(okT_exact), sum(okT_round));
        
        % 3) Also check CUES vs cue_epoch_sec 
        cue_left_sec  = cue_epoch_sec(1);
        cue_right_sec = cue_epoch_sec(2);
        left_samp_c   = cue_left_sec  * sr;   % negative
        right_samp_c  = cue_right_sec * sr;   % positive
        
        isCue   = ismember(keys, toKey(CUE_ALL));
        cueIdx  = find(isCue);
        evLatExactC = double([EEG.event(cueIdx).latency]);
        okC_exact   = (evLatExactC + left_samp_c >= 1) & (evLatExactC + right_samp_c <= EEG.pnts);
        
        % If no target is feasible, stop with a useful hint
        if ~any(okT_exact)
            earliestTarget_s = (min(evLatExactT)-1)/sr;
            fprintf(['[Hint] Earliest target at %.3f s from file start. With left=%.3f s, ' ...
                     'you need targets after %.3f s. Consider relaxing left bound.\n'], ...
                    earliestTarget_s, targ_left_sec, -targ_left_sec);
            error('No targets survive the current epoch limits.');
        end
        
        % save continues EEG in EEG_cont to use later for Master file 
        EEG_cont = EEG;  

        % --- Ensure subject + output subfolder exist BEFORE any CSV writing
        if ~exist('outEpochs','var') || isempty(outEpochs), outEpochs = pwd; end
        keepDir = fullfile(outEpochs, 'kept_indices'); 
        if ~exist(keepDir,'dir'), mkdir(keepDir); end
        
        subj = EEG_cont.subject; 
        if isempty(subj), subj = EEG_cont.setname; end
        subject = string(matlab.lang.makeValidName(string(subj)));


        %% ---------------- Target-locked epochs (padded) + labels + baselines --------
        
        % EEGLAB wants a cellstr of event types; make a safe converter
        asCellStr = @(x) cellstr(string(x(:)));
        
        EEG_targ = pop_epoch(EEG, asCellStr(TARG_ALL), targ_epoch_sec, 'epochinfo','yes');
        
        m = EEG_targ.trials;
        cue_side  = strings(m,1);  cue_soa  = strings(m,1);
        targ_side = strings(m,1);  targ_soa = strings(m,1);
        validity  = strings(m,1);
        
        bl_precue = nan(m,2);      % per-epoch pre-cue baseline mapped into target time
        
        ana_t1_ms = zeros(m,1);          % analyze from 0 ms ...
        ana_t2_ms = 800 * ones(m,1);     % ... to +800 ms
        
        for k = 1:m
            % Find the target-at-zero event index for this epoch
            evIdx = EEG_targ.epoch(k).event;
            lat   = EEG_targ.epoch(k).eventlatency;
            if iscell(lat),   lat   = cellfun(@(x) x, lat);   end
            if iscell(evIdx), evIdx = cellfun(@(x) x, evIdx); end
            z = find(abs(lat) < 1e-6, 1); if isempty(z), continue; end
            evIdxEpoch = evIdx(z);                          % index into EEG_targ.event
            if iscell(evIdxEpoch), evIdxEpoch = evIdxEpoch{1}; end
        
            % Map back to the trial table via UREVENT if available
            ii = [];
            if isfield(EEG_targ.event,'urevent') && ~isempty(EEG_targ.event(evIdxEpoch).urevent) ...
                    && isfield(TT,'targ_urev')
                u0 = EEG_targ.event(evIdxEpoch).urevent;
                ii = find([TT.targ_urev] == u0, 1);
            end
        
            if isempty(ii)
                % Fallback labels directly from the target code (robust to UNCUED)
                typ = EEG_targ.event(evIdxEpoch).type;
                cue_side(k)  = "N"; cue_soa(k) = "";
                if on(typ, TARG_L_ALL), targ_side(k) = "L";
                elseif on(typ, TARG_R_ALL), targ_side(k) = "R";
                else, targ_side(k) = "N"; end
        
                if on(typ, [TARG_L_valid_short TARG_L_invalid_short TARG_R_valid_short TARG_R_invalid_short])
                    targ_soa(k) = "short";
                elseif on(typ, [TARG_L_valid_long TARG_L_invalid_long TARG_R_valid_long TARG_R_invalid_long])
                    targ_soa(k) = "long";
                else
                    targ_soa(k) = "";
                end
        
                if on(typ, [TARG_L_valid_short TARG_L_valid_long TARG_R_valid_short TARG_R_valid_long])
                    validity(k) = "valid";
                elseif on(typ, [TARG_L_invalid_short TARG_L_invalid_long TARG_R_invalid_short TARG_R_invalid_long])
                    validity(k) = "invalid";
                else
                    validity(k) = "no_cue";
                end
                continue
            end
        
            row = TT(ii);
            cue_side(k)=row.cue_side;   cue_soa(k)=row.cue_soa;
            targ_side(k)=row.targ_side;  % keep mapping identical
            targ_soa(k)=row.targ_soa;
            validity(k)=row.validity;
        
            % SAME physiology baseline mapped into target time: [-SOA-500, -SOA-200]
            d_pts  = row.targ_lat - row.cue_lat;
            SOA_ms = 1000 * (d_pts / sr);
            bl_precue(k,:) = [-(SOA_ms + 500),  -(SOA_ms + 200)];
        
        end
        
        % ===== Patch TF baseline for UNCUED target trials to match long-SOA offset =====
        isLong = arrayfun(@(x) ~isempty(x.targ_ev) && strcmpi(x.targ_soa,'long'), TT);
        if any(isLong)
            d_pts_long  = [TT(isLong).targ_lat] - [TT(isLong).cue_lat];
            SOA_long_ms = median(1000 * (d_pts_long / sr));   % robust
        else
            SOA_long_ms = 250;
        end
        long_base = [-(SOA_long_ms + 500),  -(SOA_long_ms + 200)]; % e.g., [-750, -450] ms
        
        isNoCueEpoch = strcmpi(validity, 'no_cue') | strcmpi(cue_side, 'N');
        idxNoCue = find(isNoCueEpoch);
        
        left_ms  = EEG_targ.xmin * 1000;
        right_ms = EEG_targ.xmax * 1000;
        if ~isempty(idxNoCue) && (long_base(1) < left_ms || long_base(2) > right_ms)
            warning('Uncued baseline [%g %g] ms is outside epoch span [%g %g] ms. Consider extending targ_epoch_sec.', ...
                    long_base(1), long_base(2), left_ms, right_ms);
        end
        
        % ---------- Attach metadata to TARGET-locked epochs  ----------
        m = EEG_targ.trials;
        cue_side_c   = cellstr(cue_side);
        cue_soa_c    = cellstr(cue_soa);
        targ_side_c  = cellstr(targ_side);
        targ_soa_c   = cellstr(targ_soa);
        validity_c   = cellstr(validity);
        blp1_c       = num2cell(bl_precue(:,1));
        blp2_c       = num2cell(bl_precue(:,2));
        
        if m > 0
            [EEG_targ.epoch(1:m).cue_side]     = deal(cue_side_c{:});
            [EEG_targ.epoch(1:m).cue_soa]      = deal(cue_soa_c{:});
            [EEG_targ.epoch(1:m).targ_side]    = deal(targ_side_c{:});
            [EEG_targ.epoch(1:m).targ_soa]     = deal(targ_soa_c{:});
            [EEG_targ.epoch(1:m).validity]     = deal(validity_c{:});
            [EEG_targ.epoch(1:m).bl_precue_t1] = deal(blp1_c{:});
            [EEG_targ.epoch(1:m).bl_precue_t2] = deal(blp2_c{:});

        end
        
        for k = idxNoCue(:).'
            EEG_targ.epoch(k).bl_precue_t1 = long_base(1);
            EEG_targ.epoch(k).bl_precue_t2 = long_base(2);
        end
        fprintf('Patched uncued target epochs with long-SOA baseline [%d %d] ms (n=%d).\n', ...
                round(long_base(1)), round(long_base(2)), numel(idxNoCue));
        
        %% ---------------- Cue-locked epochs (CATCH ONLY) + labels + guard mask ----------
        % Use only catch cues so cue-locked alpha isn't contaminated by targets
        CUE_FOR_EPOCHS = [CUE_L_CATCH CUE_R_CATCH];
        
        % epoch 
        EEG_cue = pop_epoch(EEG, cellstr(string(CUE_FOR_EPOCHS(:))), cue_epoch_sec, 'epochinfo','yes');
        
        if EEG_cue.trials==0
            warning('[%s] No catch-cue epochs; skipping cue-locked block.', files(i).name);
        else
            m = EEG_cue.trials;
            cue_side  = strings(m,1); cue_soa = strings(m,1);
            targ_side = strings(m,1); targ_soa = strings(m,1);
            validity  = strings(m,1);
            targ_present = false(m,1);
        
            ana_t1_ms = zeros(m,1);
            ana_t2_ms = 800 * ones(m,1);  % with catch trials there's no target, stays 0–800 ms
        
            for k = 1:m
                evIdx = EEG_cue.epoch(k).event;
                lat   = EEG_cue.epoch(k).eventlatency;
                if iscell(lat),   lat   = cellfun(@(x) x, lat);   end
                if iscell(evIdx), evIdx = cellfun(@(x) x, evIdx); end
                z = find(abs(lat) < 1e-6, 1); if isempty(z), continue; end
                evIdxEpoch = evIdx(z); if iscell(evIdxEpoch), evIdxEpoch = evIdxEpoch{1}; end
        
                % Map via cue UREVENT if available
                ii = [];
                if isfield(EEG_cue.event,'urevent') && ~isempty(EEG_cue.event(evIdxEpoch).urevent) ...
                        && isfield(TT,'cue_urev')
                    u0 = EEG_cue.event(evIdxEpoch).urevent;
                    ii = find([TT.cue_urev] == u0, 1);
                end
                if isempty(ii), continue; end
        
                row = TT(ii);
                cue_side(k)=row.cue_side;   cue_soa(k)=row.cue_soa;   % should be "catch"
                targ_present(k)=~isempty(row.targ_ev);                 % expected false on catch
                targ_side(k)=row.targ_side; 
                targ_soa(k)=row.targ_soa;
                validity(k)=row.validity;   % "cue_only" for catch
        
                % If a target existed (shouldn't for catch), guard would shorten ana_t2_ms
                if ~isempty(row.targ_ev)
                    d_pts = row.targ_lat - row.cue_lat;
                    d_ms  = 1000 * (d_pts / sr);
                    ana_t2_ms(k) = max(0, min(800, d_ms - guard_ms));
                end
            end
        
            % Attach metadata to CUE-locked epochs
            cue_side_c     = reshape(cellstr(cue_side), [], 1);
            cue_soa_c      = reshape(cellstr(cue_soa), [], 1);
            targ_present_c = reshape(num2cell(logical(targ_present)), [], 1);
            targ_side_c    = reshape(cellstr(targ_side), [], 1);
            targ_soa_c     = reshape(cellstr(targ_soa), [], 1);
            validity_c     = reshape(cellstr(validity), [], 1);
            ana_t1_c       = reshape(num2cell(ana_t1_ms), [], 1);
            ana_t2_c       = reshape(num2cell(ana_t2_ms), [], 1);
            blc1_c         = reshape(num2cell(-500 * ones(m,1)), [], 1);
            blc2_c         = reshape(num2cell(-200 * ones(m,1)), [], 1);
        
            [EEG_cue.epoch(1:m).cue_side]     = deal(cue_side_c{:});
            [EEG_cue.epoch(1:m).cue_soa]      = deal(cue_soa_c{:});
            [EEG_cue.epoch(1:m).targ_present] = deal(targ_present_c{:});
            [EEG_cue.epoch(1:m).targ_side]    = deal(targ_side_c{:});
            [EEG_cue.epoch(1:m).targ_soa]     = deal(targ_soa_c{:});
            [EEG_cue.epoch(1:m).validity]     = deal(validity_c{:});
            [EEG_cue.epoch(1:m).ana_t1_ms]    = deal(ana_t1_c{:});
            [EEG_cue.epoch(1:m).ana_t2_ms]    = deal(ana_t2_c{:});
            [EEG_cue.epoch(1:m).bl_precue_t1] = deal(blc1_c{:});
            [EEG_cue.epoch(1:m).bl_precue_t2] = deal(blc2_c{:});
        end
 


        %% ===== Trial rejection (cue first, then target) =====
        doTrialGUI        = true;
        whichROI_cue      = "all";        % "posterior" | "parietal" | "frontal" | "all"
        whichROI_targ     = "all";
        useShadedPremark  = true;         % show autos as light-red rectangles that can be deleted
        
        % --- ROIs ---
        ROI_posterior = {'O1','Oz','O2','PO7','PO8','PO3','PO4'};
        ROI_parietal  = {'P1','P2','P3','P4','POz','P5','P6','P7','P8'};
        ROI_frontal   = {'F1','F2','F3','F4','F5','F6','F7','F8',...
                         'FC1','FC2','FC3','FC4','FC5','FC6','AF3','AF4','AF7','AF8','Fp1','Fp2'};
        
        if strcmpi(whichROI_cue,'posterior'), roi_cue = ROI_posterior;
        elseif strcmpi(whichROI_cue,'parietal'), roi_cue = ROI_parietal;
        elseif strcmpi(whichROI_cue,'frontal'),  roi_cue = ROI_frontal;
        else, roi_cue = unique([ROI_posterior ROI_parietal ROI_frontal]); end
        
        if strcmpi(whichROI_targ,'posterior'), roi_targ = ROI_posterior;
        elseif strcmpi(whichROI_targ,'parietal'), roi_targ = ROI_parietal;
        elseif strcmpi(whichROI_targ,'frontal'),  roi_targ = ROI_frontal;
        else, roi_targ = unique([ROI_posterior ROI_parietal ROI_frontal]); end
        
        % --- thresholds (lenient) ---
        TH.amp_ptp_abs_uV = 200;
        TH.amp_ptp_z      = 6.0;
        TH.jump_uV        = 80;
        TH.flat_std_uV    = 0.5;
        TH.hf_band_Hz     = [40 80];
        TH.lf_band_Hz     = [ 8 20];
        TH.hf_ratio_abs   = 3.0;
        TH.hf_ratio_z     = 3.5;
        
        set(0,'DefaultFigureWindowStyle','normal');
        set(0,'DefaultFigureVisible','on');
        
        %% ============================================================
        %% ================== CUE-locked (catch-only) =================
        %% ============================================================
        EEG = EEG_cue; fs = EEG.srate; times = EEG.times; m = EEG.trials; pnts = EEG.pnts;
        allLabs = upper(string({EEG.chanlocs.labels}));
        rix = find(ismember(allLabs, upper(string(roi_cue))));
        if isempty(rix), rix = 1:EEG.nbchan; end
        
        % per-epoch analysis window 0..ana_t2_ms (fallback 0..min(800,max(times)))
        ana_t1 = zeros(m,1);
        if isfield(EEG.epoch,'ana_t1_ms'), ana_t1 = double([EEG.epoch.ana_t1_ms]).'; end
        ana_t2 = 800*ones(m,1);
        if isfield(EEG.epoch,'ana_t2_ms'), ana_t2 = double([EEG.epoch.ana_t2_ms]).'; end
        bmask_fixed = (times>=-500 & times<=-200);
        
        PTP = nan(m,numel(rix)); JUMP = nan(m,numel(rix)); FLAT = nan(m,numel(rix)); HFR = nan(m,1);
        for k = 1:m
            tmask = (times>=ana_t1(k) & times<=ana_t2(k));
            if ~any(tmask), tmask = (times>=0 & times<=min(800, max(times))); end
        
            for c = 1:numel(rix)
                x = double(squeeze(EEG.data(rix(c),:,k)));
                if any(bmask_fixed), PTP(k,c) = max(x(bmask_fixed)) - min(x(bmask_fixed)); end
                if nnz(tmask)>=2
                    xa = x(tmask);
                    JUMP(k,c) = max(abs(diff(xa)));
                    FLAT(k,c) = std(xa);
                end
            end
        
            rr = nan(1,numel(rix));
            for c = 1:numel(rix)
                x = double(squeeze(EEG.data(rix(c),tmask,k)));
                if numel(x) > fs/5
                    [pxx,f] = pwelch(x,[],[],[],fs);
                    ph = sum(pxx(f>=TH.hf_band_Hz(1) & f<=TH.hf_band_Hz(2)));
                    pl = sum(pxx(f>=TH.lf_band_Hz(1) & f<=TH.lf_band_Hz(2)));
                    rr(c) = ph/max(pl,eps);
                end
            end
            HFR(k) = median(rr,'omitnan');
        end
        
        % robust z for PTP (per-channel) and HFR
        zPTP = nan(size(PTP));
        for c = 1:size(PTP,2)
            x = PTP(:,c); medx = median(x,'omitnan'); madx = median(abs(x-medx),'omitnan');
            zPTP(:,c) = (x - medx) ./ max(1e-9, 1.4826*madx);
        end
        x = HFR; medx = median(x,'omitnan'); madx = median(abs(x-medx),'omitnan');
        zHFR = (x - medx) ./ max(1e-9, 1.4826*madx);
        
        bad_ptp  = (PTP > TH.amp_ptp_abs_uV) | (zPTP > TH.amp_ptp_z);
        bad_jump = (JUMP > TH.jump_uV);
        bad_flat = (FLAT < TH.flat_std_uV);
        bad_hfr  = (HFR > TH.hf_ratio_abs) & (zHFR > TH.hf_ratio_z);
        
        rev = false(m,1); drop = false(m,1);
        for k = 1:m
            nBad = sum(bad_ptp(k,:) | bad_jump(k,:) | bad_flat(k,:), 'omitnan');
            if nBad >= max(2, ceil(numel(rix)/2)) || any(bad_hfr(k))
                drop(k) = true;
            elseif nBad > 0
                rev(k) = true;
            end
        end
        premark_cue = (rev | drop);
        
        % events and optional shaded rectangles
        ev = struct('latency',{},'type',{});
        for k = 1:m
            offset = (k-1)*pnts;
            t2s = @(tms) max(1, min(pnts, 1 + round(((tms/1000) - EEG.xmin) * fs)));
            ev(end+1).latency = offset + t2s(-500); ev(end).type = 'BASE_START'; %#ok<SAGROW>
            ev(end+1).latency = offset + t2s(-200); ev(end).type = 'BASE_END';
            ev(end+1).latency = offset + t2s(ana_t2(k)); ev(end).type = 'ANA_END';
            if premark_cue(k), ev(end+1).latency = offset + t2s((EEG.xmin*1000+EEG.xmax*1000)/2); ev(end).type = 'SUGG'; end
        end
  

        % ===== CUE GUI + selection -> epoch mask =====
        % check that trial2eegplot is in path
        if exist('trial2eegplot','file') ~= 2
            warning('trial2eegplot not found; shaded premarks disabled.');
            winrej_auto_cue = [];   % for the CUE block
        else
            % existing call to trial2eegplot(...)
        end
        
        % prepare the rejection mask
        if useShadedPremark
            % was: premark(:).'
            rej  = premark_cue(:).';                   % 1×m logical for CUE
            rejE = repmat(rej, EEG_cue.nbchan, 1);
            winrej_auto_cue = trial2eegplot(rej, rejE, EEG_cue.pnts, [1 0.85 0.85]);
        else
            winrej_auto_cue = [];
        end
        
        plotTitleCue = sprintf('%s | Mark bad CUE epochs', EEG_cue.setname);
        if evalin('base','exist(''CUE_WINREJ_FINAL'',''var'')'), evalin('base','clear CUE_WINREJ_FINAL'); end
        
        EEG = EEG_cue;
        eegplot(EEG.data, ...
            'srate',     EEG.srate, ...
            'eloc_file', EEG.chanlocs, ...
            'events',    ev, ...                             
            'limits',    [EEG.xmin*1000 EEG.xmax*1000], ...
            'winlength', (EEG.xmax - EEG.xmin), ...
            'dispchans', EEG.nbchan, ...
            'spacing',   50, ...
            'trialstag', 1, ...
            'winrej',    winrej_auto_cue, ...
            'butlabel',  'Save selection & close', ...
            'title',     plotTitleCue, ...
            'command',   'assignin(''base'',''CUE_WINREJ_FINAL'', TMPREJ); close(gcbf);');
        
        
        % --- Wait for window to close
        drawnow; pause(0.1);
        hFig = findall(0,'type','figure','tag','EEGPLOT','Name',plotTitleCue);
        if isempty(hFig)
            figs = findall(0,'type','figure','tag','EEGPLOT');
            if ~isempty(figs)
                hFig = figs(1);
            end
        end
        if ~isempty(hFig)
            waitfor(hFig);
        end
        
        %% --- Retrieve final rectangles from GUI
        savedViaButton = evalin('base','exist(''CUE_WINREJ_FINAL'',''var'')') == 1;
        
        winrej_final = [];
        if savedViaButton
            winrej_final = evalin('base','CUE_WINREJ_FINAL');
            evalin('base','clear CUE_WINREJ_FINAL');
        else
            % User closed without Save → treat as "no edits".
            fallbackAutosWhenNoSave = true;
            if fallbackAutosWhenNoSave
                winrej_final = winrej_auto_cue;
            else
                winrej_final = [];
            end
        end
        
        % Keep only [start end], round to samples
        if ~isempty(winrej_final)
            winrej_final = round(winrej_final(:,1:2));
        end
        
        % Prepare autos for fair "same/no edits" comparison
        auto2  = []; if ~isempty(winrej_auto_cue), auto2  = round(winrej_auto_cue(:,1:2)); end
        final2 = []; if ~isempty(winrej_final),    final2 = round(winrej_final(:,1:2));    end
        
        sameNoEdits = (size(auto2,1)==size(final2,1)) && isequal(sortrows(auto2), sortrows(final2));
        
        acceptAutosWhenNoEdits        = true;   % apply autos if window closed or identical
        acceptAutosWhenUserClearsAll  = false;  % if user saved EMPTY → keep all
        
        m_cue = EEG_cue.trials;
        pnts  = EEG_cue.pnts;
        nTot  = m_cue * pnts;
        
        fprintf('[CUE] autos=%d  GUI_rectangles=%d  (GUI edits=%d)\n', ...
            size(auto2,1), size(final2,1), ~sameNoEdits);
        
        % --- Build epoch rejection mask using overlap-per-epoch with a minimum overlap guard
        rejEpoch_cue = false(1, m_cue);
        
        if savedViaButton
            if isempty(final2)
                if acceptAutosWhenUserClearsAll
                    rejEpoch_cue = premark_cue;
                else
                    rejEpoch_cue(:) = false;
                end
            elseif sameNoEdits && acceptAutosWhenNoEdits
                rejEpoch_cue = premark_cue;
            else
                % Clean & clamp intervals
                W = final2;
                good = all(isfinite(W),2) & (W(:,2) >= W(:,1));
                W = W(good,:);
                W(:,1) = max(1, min(nTot, W(:,1)));
                W(:,2) = max(1, min(nTot, W(:,2)));
        
                % Minimum overlap threshold (samples) to count an epoch as rejected
                min_overlap = 10;  % ~20 ms at 500 Hz
        
                % For each rectangle, map to epochs via actual overlap length
                for r = 1:size(W,1)
                    s1 = W(r,1); s2 = W(r,2);
        
                    % candidate epoch index range the rectangle can touch
                    e1 = floor((s1-1)/pnts) + 1;
                    e2 = floor((s2-1)/pnts) + 1;
                    e1 = max(1, min(m_cue, e1));
                    e2 = max(1, min(m_cue, e2));
        
                    for e = e1:e2
                        ep_start = (e-1)*pnts + 1;
                        ep_end   = e*pnts;
                        overlap  = max(0, min(s2, ep_end) - max(s1, ep_start) + 1);
                        if overlap >= min_overlap
                            rejEpoch_cue(e) = true;
                        end
                    end
                end
            end
        else
            % Closed without save → "no edits"
            if acceptAutosWhenNoEdits
                rejEpoch_cue = premark_cue;
            else
                rejEpoch_cue(:) = false;
            end
        end

        
        EEG_cue.reject.rejmanual  = rejEpoch_cue;
        EEG_cue.reject.rejmanualE = repmat(rejEpoch_cue, EEG_cue.nbchan, 1);
        fprintf('[CUE] Rejected epochs: %d / %d\n', sum(rejEpoch_cue), m_cue);
        
        
        
        %% ============================================================
        %% ================= TARGET-locked  =================
        %% ============================================================
        EEG = EEG_targ; fs = EEG.srate; times = EEG.times; m = EEG.trials; pnts = EEG.pnts;
        allLabs = upper(string({EEG.chanlocs.labels}));
        rix = find(ismember(allLabs, upper(string(roi_targ))));
        if isempty(rix), rix = 1:EEG.nbchan; end
        
        tmask_const = (times>=0 & times<=800);
        hasBL = isfield(EEG.epoch,'bl_precue_t1') & isfield(EEG.epoch,'bl_precue_t2'); 
        
        PTP = nan(m,numel(rix)); JUMP = nan(m,numel(rix)); FLAT = nan(m,numel(rix)); HFR = nan(m,1);
        for k = 1:m
            if hasBL
                b1 = double(EEG.epoch(k).bl_precue_t1);
                b2 = double(EEG.epoch(k).bl_precue_t2);
                bmask = (times>=b1 & times<=b2);
            else
                bmask = (times>=-500 & times<=-200);
            end
        
            for c = 1:numel(rix)
                x = double(squeeze(EEG.data(rix(c),:,k)));
                if any(bmask), PTP(k,c) = max(x(bmask)) - min(x(bmask)); end
                xa = x(tmask_const);
                if numel(xa) >= 2
                    JUMP(k,c) = max(abs(diff(xa)));
                    FLAT(k,c) = std(xa);
                end
            end
        
            rr = nan(1,numel(rix));
            for c = 1:numel(rix)
                x = double(squeeze(EEG.data(rix(c), tmask_const, k)));
                if numel(x) > fs/5
                    [pxx,f] = pwelch(x,[],[],[],fs);
                    ph = sum(pxx(f>=TH.hf_band_Hz(1) & f<=TH.hf_band_Hz(2)));
                    pl = sum(pxx(f>=TH.lf_band_Hz(1) & f<=TH.lf_band_Hz(2)));
                    rr(c) = ph/max(pl,eps);
                end
            end
            HFR(k) = median(rr,'omitnan');
        end
        
        % robust z
        zPTP = nan(size(PTP));
        for c = 1:size(PTP,2)
            x = PTP(:,c); medx = median(x,'omitnan'); madx = median(abs(x-medx),'omitnan');
            zPTP(:,c) = (x - medx) ./ max(1e-9, 1.4826*madx);
        end
        x = HFR; medx = median(x,'omitnan'); madx = median(abs(x-medx),'omitnan');
        zHFR = (x - medx) ./ max(1e-9, 1.4826*madx);
        
        bad_ptp  = (PTP > TH.amp_ptp_abs_uV) | (zPTP > TH.amp_ptp_z);
        bad_jump = (JUMP > TH.jump_uV);
        bad_flat = (FLAT < TH.flat_std_uV);
        bad_hfr  = (HFR > TH.hf_ratio_abs) & (zHFR > TH.hf_ratio_z);
        
        rev = false(m,1); drop = false(m,1);
        for k = 1:m
            nBad = sum(bad_ptp(k,:) | bad_jump(k,:) | bad_flat(k,:), 'omitnan');
            if nBad >= max(2, ceil(numel(rix)/2)) || any(bad_hfr(k))
                drop(k) = true;
            elseif nBad > 0
                rev(k) = true;
            end
        end
        premark_targ = (rev | drop);
        
        % events (baseline per epoch) + ZERO/ANA_END + SUGG at 400 ms
        ev = struct('latency',{},'type',{});
        for k = 1:m
            offset = (k-1)*pnts;
            t2s = @(tms) max(1, min(pnts, 1 + round(((tms/1000) - EEG.xmin) * fs)));
            if hasBL
                b1 = double(EEG.epoch(k).bl_precue_t1);
                b2 = double(EEG.epoch(k).bl_precue_t2);
            else
                b1 = -500; b2 = -200;
            end
            ev(end+1).latency = offset + t2s(b1);   ev(end).type = 'BASE_START'; %#ok<SAGROW>
            ev(end+1).latency = offset + t2s(b2);   ev(end).type = 'BASE_END';
            ev(end+1).latency = offset + t2s(0);    ev(end).type = 'ZERO';
            ev(end+1).latency = offset + t2s(800);  ev(end).type = 'ANA_END';
            if premark_targ(k), ev(end+1).latency = offset + t2s(400); ev(end).type = 'SUGG'; end
        end
        
        % ===== TARGET GUI + selection -> epoch mask =====
        % check that trial2eegplot is on path
        if exist('trial2eegplot','file') ~= 2
            warning('trial2eegplot not found; shaded premarks disabled.');
            winrej_auto_targ = [];   % for the target block
        else
            % existing call to trial2eegplot(...)
        end


        if useShadedPremark
            % was: premark(:).'
            rej  = premark_targ(:).';                  % 1×m logical for TARGET
            rejE = repmat(rej, EEG_targ.nbchan, 1);
            winrej_auto_targ = trial2eegplot(rej, rejE, EEG_targ.pnts, [1 0.85 0.85]);
        else
            winrej_auto_targ = [];
        end
        
        plotTitleTarg = sprintf('%s | Mark bad TARGET epochs', EEG_targ.setname);
        if evalin('base','exist(''TARG_WINREJ_FINAL'',''var'')'), evalin('base','clear TARG_WINREJ_FINAL'); end
        
        EEG = EEG_targ;
        eegplot(EEG.data, ...
            'srate',     EEG.srate, ...
            'eloc_file', EEG.chanlocs, ...
            'events',    ev, ...                              
            'limits',    [EEG.xmin*1000 EEG.xmax*1000], ...
            'winlength', (EEG.xmax - EEG.xmin), ...
            'dispchans', EEG.nbchan, ...
            'spacing',   50, ...
            'trialstag', 1, ...
            'winrej',    winrej_auto_targ, ...
            'butlabel',  'Save selection & close', ...
            'title',     plotTitleTarg, ...
            'command',   'assignin(''base'',''TARG_WINREJ_FINAL'', TMPREJ); close(gcbf);');
        

        fprintf('[TARGET] autos=%d  rectangles=%d\n', sum(premark_targ), size(winrej_auto_targ,1));
        
        % --- Wait for window to close
        drawnow; pause(0.1);
        hFig = findall(0,'type','figure','tag','EEGPLOT','Name',plotTitleTarg);
        if isempty(hFig)
            figs = findall(0,'type','figure','tag','EEGPLOT');
            if ~isempty(figs)
                hFig = figs(1);
            end
        end
        if ~isempty(hFig)
            waitfor(hFig);
        end
        
        %% --- Retrieve final rectangles from GUI 
        savedViaButton = evalin('base','exist(''TARG_WINREJ_FINAL'',''var'')') == 1;
        
        winrej_final = [];
        if savedViaButton
            winrej_final = evalin('base','TARG_WINREJ_FINAL');
            evalin('base','clear TARG_WINREJ_FINAL');
        else
            % User closed without pressing Save -> treat as "no edits".
            fallbackAutosWhenNoSave = true;
            if fallbackAutosWhenNoSave
                winrej_final = winrej_auto_targ;
            else
                winrej_final = [];
            end
        end
        
        % Keep only [start end], round to samples
        if ~isempty(winrej_final)
            winrej_final = round(winrej_final(:,1:2));
        end
        
        % Prepare autos & finals for fair "no edits" comparison (order-insensitive)
        auto2  = []; if ~isempty(winrej_auto_targ), auto2  = round(winrej_auto_targ(:,1:2)); end
        final2 = []; if ~isempty(winrej_final),     final2 = round(winrej_final(:,1:2));     end
        sameNoEdits = (size(auto2,1)==size(final2,1)) && isequal(sortrows(auto2), sortrows(final2));
        
        acceptAutosWhenNoEdits       = true;   % apply autos if window closed or identical
        acceptAutosWhenUserClearsAll = false;  % if user saved EMPTY → keep all
        
        % report 
        fprintf('[TARGET] autos=%d  GUI_rectangles=%d  (GUI edits=%d)\n', ...
            size(auto2,1), size(final2,1), ~sameNoEdits);
        
        % --- Build epoch rejection mask using overlap-per-epoch with a minimum overlap guard
        m_targ     = EEG_targ.trials;
        pnts_targ  = EEG_targ.pnts;
        nTot_targ  = m_targ * pnts_targ;
        rejEpoch_targ = false(1, m_targ);
        
        if savedViaButton
            if isempty(final2)
                % User explicitly saved empty set -> keep all 
                if acceptAutosWhenUserClearsAll
                    rejEpoch_targ = premark_targ;
                else
                    rejEpoch_targ(:) = false;
                end
            elseif sameNoEdits && acceptAutosWhenNoEdits
                % Same as autos -> trust autos
                rejEpoch_targ = premark_targ;
            else
                % User edited -> map rectangles to epochs with overlap threshold
                W = final2;
                good = all(isfinite(W),2) & (W(:,2) >= W(:,1));
                W = W(good,:);
                if ~isempty(W)
                    % clamp to valid sample range
                    W(:,1) = max(1, min(nTot_targ, W(:,1)));
                    W(:,2) = max(1, min(nTot_targ, W(:,2)));
        
                    % Minimum required overlap (samples) for an epoch to count as rejected
                    min_overlap = 10;   % ≈20 ms at 500 Hz; adjust as you like
        
                    for r = 1:size(W,1)
                        s1 = W(r,1); s2 = W(r,2);
        
                        % candidate epoch indices touched by this rectangle
                        e1 = floor((s1-1)/pnts_targ) + 1;
                        e2 = floor((s2-1)/pnts_targ) + 1;
                        e1 = max(1, min(m_targ, e1));
                        e2 = max(1, min(m_targ, e2));
        
                        for e = e1:e2
                            ep_start = (e-1)*pnts_targ + 1;
                            ep_end   = e*pnts_targ;
                            overlap  = max(0, min(s2, ep_end) - max(s1, ep_start) + 1);
                            if overlap >= min_overlap
                                rejEpoch_targ(e) = true;
                            end
                        end
                    end
                end
            end
        else
            % Closed without Save -> "no edits" policy
            if acceptAutosWhenNoEdits
                rejEpoch_targ = premark_targ;
            else
                rejEpoch_targ(:) = false;
            end
        end

        
        EEG_targ.reject.rejmanual  = rejEpoch_targ;
        EEG_targ.reject.rejmanualE = repmat(rejEpoch_targ, EEG_targ.nbchan, 1);
        
        fprintf('[TARGET] Rejected epochs: %d / %d\n', sum(rejEpoch_targ), m_targ);





       % ===================== QA MASTER (UNIFIED: PREFLIGHT + POST-REJ) =====================
        
        if ~exist('outEpochs','var') || isempty(outEpochs), outEpochs = pwd; end
        qaCsvPath = fullfile(outEpochs, 'posner_QA_MASTER.csv');
        dictPath  = fullfile(outEpochs, 'posner_QA_DICTIONARY.csv');
        keepDir   = fullfile(outEpochs, 'kept_indices'); if ~exist(keepDir,'dir'), mkdir(keepDir); end
        
        % ---------- Subject name ----------
        subj = EEG_cont.subject; if isempty(subj), subj = EEG_cont.setname; end
        subject = string(matlab.lang.makeValidName(string(subj)));
        
        % ---------- Helpers ----------
        toKey = @(x) regexprep(upper(string(x)),'\D','');   % keep only digits
        asCellStr = @(x) cellstr(string(x(:)));
        
        % ---------- PREFLIGHT (recompute robustly if prior vars not in scope) ----------
        % Event code groups 
        CUE_L_short = {'110'}; CUE_R_short = {'120'};
        CUE_L_long  = {'210'}; CUE_R_long  = {'220'};
        CUE_L_CATCH = {'10'};  CUE_R_CATCH = {'20'};
        CUE_ALL = [CUE_L_short CUE_R_short CUE_L_long CUE_R_long CUE_L_CATCH CUE_R_CATCH];
        
        TARG_L_valid_long    = {'211'};  TARG_L_valid_short   = {'111'};
        TARG_L_invalid_long  = {'221'};  TARG_L_invalid_short = {'121'};
        TARG_R_valid_long    = {'222'};  TARG_R_valid_short   = {'122'};
        TARG_R_invalid_long  = {'212'};  TARG_R_invalid_short = {'112'};
        TARG_L_uncued = {'1'}; TARG_R_uncued = {'2'};
        
        TARG_L_ALL = [TARG_L_valid_long TARG_L_valid_short TARG_L_invalid_long TARG_L_invalid_short TARG_L_uncued];
        TARG_R_ALL = [TARG_R_valid_long TARG_R_valid_short TARG_R_invalid_long TARG_R_invalid_short TARG_R_uncued];
        TARG_ALL   = [TARG_L_ALL TARG_R_ALL];
        
        % Keys on the *continuous* stream
        keys = toKey({EEG_cont.event.type});
        
        % Raw-stream counts (before epoching)
        n_cues_total     = sum(ismember(keys, toKey(CUE_ALL)));
        n_cue_L_short    = sum(ismember(keys, toKey(CUE_L_short)));
        n_cue_R_short    = sum(ismember(keys, toKey(CUE_R_short)));
        n_cue_L_long     = sum(ismember(keys, toKey(CUE_L_long)));
        n_cue_R_long     = sum(ismember(keys, toKey(CUE_R_long)));
        n_cue_L_catch    = sum(ismember(keys, toKey(CUE_L_CATCH)));
        n_cue_R_catch    = sum(ismember(keys, toKey(CUE_R_CATCH)));
        
        isTarg_all  = ismember(keys, toKey(TARG_ALL));
        isCuedTarg  = ismember(keys, toKey([TARG_L_valid_long TARG_L_valid_short TARG_L_invalid_long TARG_L_invalid_short ...
                                            TARG_R_valid_long TARG_R_valid_short TARG_R_invalid_long TARG_R_invalid_short]));
        isUncTarg   = ismember(keys, toKey([TARG_L_uncued TARG_R_uncued]));
        isValidRaw  = ismember(keys, toKey([TARG_L_valid_long TARG_L_valid_short TARG_R_valid_long TARG_R_valid_short]));
        isInvalidRaw= ismember(keys, toKey([TARG_L_invalid_long TARG_L_invalid_short TARG_R_invalid_long TARG_R_invalid_short]));
        isShortRaw  = ismember(keys, toKey([TARG_L_valid_short TARG_L_invalid_short TARG_R_valid_short TARG_R_invalid_short]));
        isLongRaw   = ismember(keys, toKey([TARG_L_valid_long  TARG_L_invalid_long  TARG_R_valid_long  TARG_R_invalid_long]));
        isLeftRaw   = ismember(keys, toKey(TARG_L_ALL));
        isRightRaw  = ismember(keys, toKey(TARG_R_ALL));
        
        n_targets_total_all    = sum(isTarg_all);
        n_targets_total_cued   = sum(isCuedTarg);
        n_targets_total_uncued = sum(isUncTarg);
        n_targets_valid_total  = sum(isValidRaw);
        n_targets_invalid_total= sum(isInvalidRaw);
        n_targets_short_total  = sum(isShortRaw);
        n_targets_long_total   = sum(isLongRaw);
        n_targets_left_total   = sum(isLeftRaw);
        n_targets_right_total  = sum(isRightRaw);
        
        % Trials total = #cues + #uncued targets
        n_trials_total = n_cues_total + n_targets_total_uncued;
        
        % Epoch windows (recompute from chosen spans)
        targ_left_sec  = EEG_targ.xmin;   targ_right_sec = EEG_targ.xmax;
        cue_left_sec   = EEG_cue.xmin;    cue_right_sec  = EEG_cue.xmax;
        targ_window_left_ms  = round(1000*abs(targ_left_sec));
        targ_window_right_ms = round(1000*targ_right_sec);
        cue_window_left_ms   = round(1000*abs(cue_left_sec));
        cue_window_right_ms  = round(1000*cue_right_sec);
        
        % Feasibility (exact) check for targets & cues w.r.t. the continuous data
        sr = EEG_cont.srate;
        targIdx = find(ismember(keys, toKey(TARG_ALL)));
        cueIdx  = find(ismember(keys, toKey(CUE_ALL)));
        evLatT  = double([EEG_cont.event(targIdx).latency]);
        evLatC  = double([EEG_cont.event(cueIdx).latency]);
        
        okT_exact = (evLatT + targ_left_sec*sr  >= 1) & (evLatT + targ_right_sec*sr <= EEG_cont.pnts);
        okC_exact = (evLatC + cue_left_sec*sr   >= 1) & (evLatC + cue_right_sec*sr  <= EEG_cont.pnts);
        
        n_targets_feasible_exact = sum(okT_exact);
        n_cues_feasible_exact    = sum(okC_exact);
        n_targets_epoched        = EEG_targ.trials;
        n_cues_epoched           = EEG_cue.trials;
        
        % Preflight drop reasons for targets (relative to chosen window)
        if isempty(targIdx)
            n_targets_drop_early=0; n_targets_drop_late=0; n_targets_drop_both=0; n_targets_drop_boundary=0;
        else
            badT  = targIdx(~okT_exact);
            needL = -targ_left_sec; needR = targ_right_sec;
            if isempty(badT)
                n_targets_drop_early=0; n_targets_drop_late=0; n_targets_drop_both=0;
            else
                haveL  = ((double([EEG_cont.event(badT).latency])-1) ./ sr).';
                haveR  = ((EEG_cont.pnts - double([EEG_cont.event(badT).latency])) ./ sr).';
                missL  = max(0, needL - haveL);
                missR  = max(0, needR - haveR);
                n_targets_drop_both  = sum(missL>0 & missR>0);
                n_targets_drop_early = sum(missL>0 & missR==0);
                n_targets_drop_late  = sum(missL==0 & missR>0);
            end
            % boundary-in-window among feasible-but-missing target epochs
            isBoundary = strcmpi({EEG_cont.event.type}, 'BOUNDARY');
            bLat = double([EEG_cont.event(isBoundary).latency]);
            feas_idx  = targIdx(okT_exact);
            % map epoched targets (urevent) to detect which feasible ones did not get epoched
            epoch_urev = nan(n_targets_epoched,1);
            if isfield(EEG_targ.event,'urevent')
                for k = 1:n_targets_epoched
                    evIdx = EEG_targ.epoch(k).event;
                    lat   = EEG_targ.epoch(k).eventlatency;
                    if iscell(lat),   lat=cellfun(@(x) x, lat);  end
                    if iscell(evIdx), evIdx=cellfun(@(x) x, evIdx); end
                    z = find(abs(lat)<1e-6,1);
                    if ~isempty(z)
                        e = evIdx(z); if iscell(e), e=e{1}; end
                        if ~(isscalar(e) && e>=1 && e<=numel(EEG.event))
                            keep(k)     = false;
                            % choose a reason field that exists in this block:
                            if exist('reason_v','var'); reason_v(k) = "bad_event_index"; end
                            continue
                        end

                        u = EEG_targ.event(e).urevent;
                        if ~isempty(u), epoch_urev(k) = double(u); end
                    end
                end
            end
            feas_urev = nan(size(feas_idx));
            if isfield(EEG_cont.event,'urevent') && ~isempty([EEG_cont.event.urevent])
                feas_urev = double([EEG_cont.event(feas_idx).urevent]);
            end
            missing_mask = ~ismember(feas_urev, epoch_urev);
            missing_idx_all = feas_idx(missing_mask);
            n_targets_drop_boundary = 0;
            if ~isempty(missing_idx_all) && ~isempty(bLat)
                left_samp_t  = targ_left_sec  * sr;
                right_samp_t = targ_right_sec * sr;
                for e = missing_idx_all(:).'
                    L = double(EEG_cont.event(e).latency);
                    winStart = L + left_samp_t;
                    winEnd   = L + right_samp_t;
                    if any( bLat >= ceil(winStart) & bLat <= floor(winEnd) )
                        n_targets_drop_boundary = n_targets_drop_boundary + 1;
                    end
                end
            end
        end
        
        % ---------- POST-REJECTION (KEPT) ----------
        % Inclusion masks BEFORE actually removing epochs
        keep_cue  = true(1, EEG_cue.trials);  if isfield(EEG_cue.reject,'rejmanual')  && ~isempty(EEG_cue.reject.rejmanual),  keep_cue  = ~EEG_cue.reject.rejmanual;  end
        keep_targ = true(1, EEG_targ.trials); if isfield(EEG_targ.reject,'rejmanual') && ~isempty(EEG_targ.reject.rejmanual), keep_targ = ~EEG_targ.reject.rejmanual; end

        
       % ---------- Build CUE kept-index table (length-safe) ----------
        m_cue = EEG_cue.trials;
        epoch_cue = (1:m_cue).';
        kept_cue  = reshape(keep_cue(:), [], 1);
        
        % Safe extractors that always return column vectors of length m_cue
        cue_side_v = strings(m_cue,1);
        if isfield(EEG_cue,'epoch') && ~isempty(EEG_cue.epoch) && isfield(EEG_cue.epoch,'cue_side')
            tmp = string({EEG_cue.epoch.cue_side}).';
            if numel(tmp) == m_cue
                cue_side_v = tmp;
            else
                warning('[CUE kept] cue_side length %d ~= trials %d; filling defaults', numel(tmp), m_cue);
            end
        end
        
        % quick length audit before table creation
        fprintf('[CUE kept] lens epoch=%d kept=%d cue_side=%d\n', numel(epoch_cue), numel(kept_cue), numel(cue_side_v));
        
        Tcue = table(epoch_cue, kept_cue, cue_side_v, ...
            'VariableNames', {'epoch','kept','cue_side'});
        writetable(Tcue, fullfile(keepDir, sprintf('%s_cue_kept_indices.csv', subject)), 'Delimiter',';');
        
        % ---------- Build TARGET kept-index table (length-safe) ----------
        m_targ = EEG_targ.trials;
        epoch_t  = (1:m_targ).';
        kept_t   = reshape(keep_targ(:), [], 1);
        
        validity_v  = strings(m_targ,1);
        targ_side_v = strings(m_targ,1);
        targ_soa_v  = strings(m_targ,1);
        
        if isfield(EEG_targ,'epoch') && ~isempty(EEG_targ.epoch)
            % validity
            if isfield(EEG_targ.epoch,'validity')
                tmp = string({EEG_targ.epoch.validity}).';
                if numel(tmp)==m_targ, validity_v = tmp;
                else, warning('[TARGET kept] validity length %d ~= trials %d; filling defaults', numel(tmp), m_targ);
                end
            end
            % targ_side
            if isfield(EEG_targ.epoch,'targ_side')
                tmp = string({EEG_targ.epoch.targ_side}).';
                if numel(tmp)==m_targ, targ_side_v = tmp;
                else, warning('[TARGET kept] targ_side length %d ~= trials %d; filling defaults', numel(tmp), m_targ);
                end
            end
            % targ_soa
            if isfield(EEG_targ.epoch,'targ_soa')
                tmp = string({EEG_targ.epoch.targ_soa}).';
                if numel(tmp)==m_targ, targ_soa_v = tmp;
                else, warning('[TARGET kept] targ_soa length %d ~= trials %d; filling defaults', numel(tmp), m_targ);
                end
            end
        end
        
        % quick length audit before table creation
        fprintf('[TARGET kept] lens epoch=%d kept=%d validity=%d side=%d soa=%d\n', ...
            numel(epoch_t), numel(kept_t), numel(validity_v), numel(targ_side_v), numel(targ_soa_v));
        
        Tt = table(epoch_t, kept_t, validity_v, targ_side_v, targ_soa_v, ...
            'VariableNames', {'epoch','kept','validity','targ_side','targ_soa'});
        writetable(Tt, fullfile(keepDir, sprintf('%s_target_kept_indices.csv', subject)), 'Delimiter',';');

        
% ========== Robust kept-counts (no gtar/gcue; align strictly to current sets) ==========
% Target side — use vectors that were constructed for Tt (length m_targ)
V  = reshape(string(validity_v),  1, []);         % "valid" | "invalid" | "no_cue"
S  = reshape(string(targ_side_v), 1, []);         % "L" | "R" | "N"
SO = reshape(string(targ_soa_v),  1, []);         % "short" | "long" | ""
K  = reshape(logical(keep_targ),  1, []);         % logical inclusion mask

% Sanity: all 1×m_targ
if ~isequal(size(V),size(S),size(SO),size(K))
    error('QA kept-counts: size mismatch (targets): V=%s S=%s SO=%s K=%s', ...
        mat2str(size(V)), mat2str(size(S)), mat2str(size(SO)), mat2str(size(K)));
end

isValid   = (V  == "valid");
isInvalid = (V  == "invalid");
isNoCue   = (V  == "no_cue");
isL       = (S  == "L");
isR       = (S  == "R");
isShort   = (SO == "short");
isLong    = (SO == "long");

n_kept_targ_valid_L_short   = sum(K & isValid   & isL & isShort);
n_kept_targ_valid_L_long    = sum(K & isValid   & isL & isLong);
n_kept_targ_valid_R_short   = sum(K & isValid   & isR & isShort);
n_kept_targ_valid_R_long    = sum(K & isValid   & isR & isLong);

n_kept_targ_invalid_L_short = sum(K & isInvalid & isL & isShort);
n_kept_targ_invalid_L_long  = sum(K & isInvalid & isL & isLong);
n_kept_targ_invalid_R_short = sum(K & isInvalid & isR & isShort);
n_kept_targ_invalid_R_long  = sum(K & isInvalid & isR & isLong);

n_kept_targ_nocue_L         = sum(K & isNoCue & isL);
n_kept_targ_nocue_R         = sum(K & isNoCue & isR);

n_kept_targ_cued_total      = sum(K & ~isNoCue);
n_kept_targ_nocue_total     = sum(K &  isNoCue);
n_kept_targ_total           = sum(K);

% Cue side — use vectors that were constructed for Tcue (length m_cue)
Sc = reshape(string(cue_side_v), 1, []);
Kc = reshape(logical(keep_cue),  1, []);

if ~isequal(size(Sc), size(Kc))
    error('QA kept-counts: size mismatch (cues): Sc=%s Kc=%s', mat2str(size(Sc)), mat2str(size(Kc)));
end

n_kept_cue_catch_L = sum(Kc & (Sc=="L"));
n_kept_cue_catch_R = sum(Kc & (Sc=="R"));
n_kept_cue_total   = sum(Kc);


        
        % ---------- Assemble one subject row (ALL columns) ----------
        Row = table( ...
            subject, ...
            n_trials_total, ...
            n_cues_total, n_cue_L_short, n_cue_R_short, n_cue_L_long, n_cue_R_long, n_cue_L_catch, n_cue_R_catch, ...
            n_targets_total_all, n_targets_total_cued, n_targets_total_uncued, ...
            n_targets_valid_total, n_targets_invalid_total, ...
            n_targets_short_total, n_targets_long_total, ...
            n_targets_left_total, n_targets_right_total, ...
            targ_window_left_ms, targ_window_right_ms, ...
            n_targets_feasible_exact, n_targets_epoched, ...
            n_targets_drop_early, n_targets_drop_late, n_targets_drop_both, n_targets_drop_boundary, ...
            cue_window_left_ms, cue_window_right_ms, ...
            n_cues_feasible_exact, n_cues_epoched, ...
            n_kept_cue_catch_L, n_kept_cue_catch_R, n_kept_cue_total, ...
            n_kept_targ_valid_L_short, n_kept_targ_valid_L_long, n_kept_targ_valid_R_short, n_kept_targ_valid_R_long, ...
            n_kept_targ_invalid_L_short, n_kept_targ_invalid_L_long, n_kept_targ_invalid_R_short, n_kept_targ_invalid_R_long, ...
            n_kept_targ_nocue_L, n_kept_targ_nocue_R, ...
            n_kept_targ_cued_total, n_kept_targ_nocue_total, n_kept_targ_total, ...
            'VariableNames', [ ...
              "subject", ...
              "n_trials_total", ...
              "n_cues_total","n_cue_L_short","n_cue_R_short","n_cue_L_long","n_cue_R_long","n_cue_L_catch","n_cue_R_catch", ...
              "n_targets_total_all","n_targets_total_cued","n_targets_total_uncued", ...
              "n_targets_valid_total","n_targets_invalid_total", ...
              "n_targets_short_total","n_targets_long_total", ...
              "n_targets_left_total","n_targets_right_total", ...
              "targ_window_left_ms","targ_window_right_ms", ...
              "n_targets_feasible_exact","n_targets_epoched", ...
              "n_targets_drop_early","n_targets_drop_late","n_targets_drop_both","n_targets_drop_boundary", ...
              "cue_window_left_ms","cue_window_right_ms", ...
              "n_cues_feasible_exact","n_cues_epoched", ...
              "n_kept_cue_catch_L","n_kept_cue_catch_R","n_kept_cue_total", ...
              "n_kept_targ_valid_L_short","n_kept_targ_valid_L_long","n_kept_targ_valid_R_short","n_kept_targ_valid_R_long", ...
              "n_kept_targ_invalid_L_short","n_kept_targ_invalid_L_long","n_kept_targ_invalid_R_short","n_kept_targ_invalid_R_long", ...
              "n_kept_targ_nocue_L","n_kept_targ_nocue_R", ...
              "n_kept_targ_cued_total","n_kept_targ_nocue_total","n_kept_targ_total" ...
            ] ...
        );
        
        % ---------- Merge/overwrite row in master CSV ----------
        if isfile(qaCsvPath)
            M = readtable(qaCsvPath, 'TextType','string', 'Delimiter',';');
        else
            M = table();
        end
        
        % Ensure all columns in Row exist in M (type-aware)
        want = Row.Properties.VariableNames;
        have = M.Properties.VariableNames;
        missing_in_M = setdiff(want, have);
        for v = string(missing_in_M)
            proto = Row.(v);
            if isstring(proto) || iscellstr(proto)
                M.(v) = strings(height(M),1);
            else
                M.(v) = zeros(height(M),1);
            end
        end
        
        % Ensure all columns in M exist in Row (type-aware)
        missing_in_Row = setdiff(have, want);
        for v = string(missing_in_Row)
            proto = M.(v);
            if isstring(proto) || iscellstr(proto)
                Row.(v) = strings(height(Row),1);
            else
                Row.(v) = zeros(height(Row),1);
            end
        end
        
        % Align column order
        M   = M(:, union(have, want, 'stable'));
        Row = Row(:, M.Properties.VariableNames);

        
        % Overwrite or append by subject
        if ~ismember("subject", string(M.Properties.VariableNames))
            M = Row;
        else
            if ~isstring(M.subject), M.subject = string(M.subject); end
            idx = find(M.subject == subject, 1);
            if isempty(idx), M = [M; Row]; else, M(idx,:) = Row; end
        end
        
        if ismember("subject", M.Properties.VariableNames)
            M = sortrows(M, 'subject');
        end
        writetable(M, qaCsvPath, 'Delimiter',';');
        fprintf('[QA] Wrote/updated unified QA master: %s\n', qaCsvPath);
        
     % ---------- Column dictionary (auto-synced; no containers.Map) ----------
dictVersion = 1;  % bump whenever texts change

Mcols = string(M.Properties.VariableNames).';   % column string array

% Known descriptions (keys and texts same length, both string arrays)
knownKeys = string([ ...
  "_version_", ...
  "subject","n_trials_total", ...
  "n_cues_total","n_cue_L_short","n_cue_R_short","n_cue_L_long","n_cue_R_long","n_cue_L_catch","n_cue_R_catch", ...
  "n_targets_total_all","n_targets_total_cued","n_targets_total_uncued", ...
  "n_targets_valid_total","n_targets_invalid_total", ...
  "n_targets_short_total","n_targets_long_total", ...
  "n_targets_left_total","n_targets_right_total", ...
  "targ_window_left_ms","targ_window_right_ms", ...
  "n_targets_feasible_exact","n_targets_epoched", ...
  "n_targets_drop_early","n_targets_drop_late","n_targets_drop_both","n_targets_drop_boundary", ...
  "cue_window_left_ms","cue_window_right_ms", ...
  "n_cues_feasible_exact","n_cues_epoched", ...
  "n_kept_cue_catch_L","n_kept_cue_catch_R","n_kept_cue_total", ...
  "n_kept_targ_valid_L_short","n_kept_targ_valid_L_long","n_kept_targ_valid_R_short","n_kept_targ_valid_R_long", ...
  "n_kept_targ_invalid_L_short","n_kept_targ_invalid_L_long","n_kept_targ_invalid_R_short","n_kept_targ_invalid_R_long", ...
  "n_kept_targ_nocue_L","n_kept_targ_nocue_R", ...
  "n_kept_targ_cued_total","n_kept_targ_nocue_total","n_kept_targ_total" ]);

knownDesc = string([ ...
  num2str(dictVersion), ...
  "Participant identifier (sanitized from EEG.subject/EEG.setname).", ...
  "Total trial *units* in raw stream = #all cues + #uncued targets.", ...
  "All cues in raw stream (any side/SOA).", ...
  "Left short cues (code 110) in raw stream.", ...
  "Right short cues (code 120) in raw stream.", ...
  "Left long cues (code 210) in raw stream.", ...
  "Right long cues (code 220) in raw stream.", ...
  "Left catch cues (code 10) in raw stream.", ...
  "Right catch cues (code 20) in raw stream.", ...
  "All targets in raw stream (cued + uncued).", ...
  "Cued targets only (valid+invalid; all sides/SOAs).", ...
  "Uncued targets only (codes 1/2).", ...
  "Valid targets in raw stream (all sides/SOAs).", ...
  "Invalid targets in raw stream (all sides/SOAs).", ...
  "Short-SOA targets in raw stream (all sides; valid+invalid).", ...
  "Long-SOA targets in raw stream (all sides; valid+invalid).", ...
  "Targets on left side in raw stream (incl. uncued left).", ...
  "Targets on right side in raw stream (incl. uncued right).", ...
  "ABS(left) extent of TARGET epoch window in ms (e.g., 1200 = −1.2 s).", ...
  "Right extent of TARGET epoch window in ms (e.g., 1300 = +1.3 s).", ...
  "Targets far enough from edges (EXACT test) given TARGET window.", ...
  "TARGET-locked epochs created by pop_epoch (before rejection).", ...
  "Feasible targets dropped preflight: too early for chosen window.", ...
  "Feasible targets dropped preflight: too late for chosen window.", ...
  "Feasible targets dropped preflight for both early & late.", ...
  "Feasible targets that failed due to a BOUNDARY in the epoch window.", ...
  "ABS(left) extent of CUE epoch window in ms.", ...
  "Right extent of CUE epoch window in ms.", ...
  "Cues far enough from edges (EXACT test) given CUE window.", ...
  "CUE-locked epochs created by pop_epoch (before rejection).", ...
  "KEPT (post-rejection) catch-cue epochs on the left.", ...
  "KEPT (post-rejection) catch-cue epochs on the right.", ...
  "KEPT (post-rejection) catch-cue epochs total.", ...
  "KEPT cued TARGETS: valid, left, short SOA.", ...
  "KEPT cued TARGETS: valid, left, long SOA.", ...
  "KEPT cued TARGETS: valid, right, short SOA.", ...
  "KEPT cued TARGETS: valid, right, long SOA.", ...
  "KEPT cued TARGETS: invalid, left, short SOA.", ...
  "KEPT cued TARGETS: invalid, left, long SOA.", ...
  "KEPT cued TARGETS: invalid, right, short SOA.", ...
  "KEPT cued TARGETS: invalid, right, long SOA.", ...
  "KEPT UNCUED targets on the left (code 1).", ...
  "KEPT UNCUED targets on the right (code 2).", ...
  "KEPT cued TARGETS total (valid+invalid).", ...
  "KEPT UNCUED targets total.", ...
  "KEPT TARGETS total (cued + uncued)." ]);

% Build description vector aligned to actual columns; fill unknowns with ""
Desc = strings(numel(Mcols),1);
[tf,loc] = ismember(Mcols, knownKeys);
Desc(tf)  = knownDesc(loc(tf));
Desc(~tf) = "";   % no text yet for these columns

% First row: version marker
Col_full  = ["_version_"; Mcols];
Desc_full = [string(num2str(dictVersion)); Desc];

Dict = table(Col_full, Desc_full, 'VariableNames', ["column","description"]);
writetable(Dict, dictPath, 'Delimiter',';');


        % save data set with marked rejection trials BEFORE dropping:
        pop_saveset(EEG_cue,  'filename',[safeSubj '_cue_marked.set'],  'filepath', outEpochs);
        pop_saveset(EEG_targ, 'filename',[safeSubj '_targ_marked.set'], 'filepath', outEpochs);

        
        % To actually drop afterwards:
        EEG_cue  = pop_rejepoch(EEG_cue,  EEG_cue.reject.rejmanual, 0);
        EEG_targ = pop_rejepoch(EEG_targ, EEG_targ.reject.rejmanual, 0);

        % Alias to *_clean for clarity
        EEG_cue_clean  = EEG_cue;
        EEG_targ_clean = EEG_targ;

        % Save analysis sets:
        pop_saveset(EEG_cue_clean,  'filename',[safeSubj '_cue_clean.set'],  'filepath', outEpochs);
        pop_saveset(EEG_targ_clean, 'filename',[safeSubj '_targ_clean.set'], 'filepath', outEpochs);

        % --- Clean up before moving to next participant
        close(findall(0,'type','figure','tag','EEGPLOT'));  % close leftover scroll windows
        clear EEG EEG_targ EEG_cue EEG_cont TT              % free large variables
        clearvars -global ALLEEG CURRENTSET EEGTMP
        eeglab redraw;    

        catch ME
        % Pick a safe label for the failing item
        if i <= numel(files) && isfield(files,'name')
            failName = files(i).name;
        else
            failName = sprintf('item_%d', i);
        end

        warning('[Stage D/E] %s failed (%s): %s', failName, ME.identifier, ME.message);


        % --- Error cleanup -----
        close(findall(0,'type','figure','tag','EEGPLOT'));  % close leftover scroll windows
        clear EEG EEG_targ EEG_cue EEG_cont TT              % free large variables
        clearvars -global ALLEEG CURRENTSET EEGTMP
        eeglab redraw;                                     % refresh EEGLAB internal state

        continue;  % be explicit about moving to the next subject
    end
end





%% 
%% ================================================================
%% 13. Behavior-only post-filter (reopen *_clean, drop FA / wrong-side)
%% ================================================================
% EEGLAB init 
addpath('...\eeglab2025.0.0');
eeglab; close;  % init quietly

inEpochs     = '...\Epochs';
outEpochsBeh = fullfile(inEpochs, 'BehaviorClean');
if ~isfolder(outEpochsBeh), mkdir(outEpochsBeh); end

set(0,'DefaultFigureWindowStyle','normal');  % keep EEGLAB plots undocked if any appear

% Collect subjects from available *_targ_clean.set files (one pass)
ft = dir(fullfile(inEpochs, '*_targ_clean.set'));
fc = dir(fullfile(inEpochs, '*_cue_clean.set'));

% Build a unique subject list from both cues and targets
names = strings(0,1);

for i = 1:numel(ft)
    [~, nm, ~] = fileparts(ft(i).name);           
    nm = regexprep(nm, '_targ_clean$', '');       
    names(end+1,1) = string(nm);
end
for i = 1:numel(fc)
    [~, nm, ~] = fileparts(fc(i).name);           
    nm = regexprep(nm, '_cue_clean$', '');        
    names(end+1,1) = string(nm);
end

% Deduplicate and drop accidental empties
subjects = unique(names(names ~= ""));
fprintf('[BehaviorFilter] Found %d subjects: %s\n', numel(subjects), strjoin(subjects.', ', '));

% Tolerance for locating the "event at 0 ms" inside an epoch (protects against FP jitter)
epsZero = 1e-6;  % tolerance in the SAME units as EEG.epoch(...).eventlatency

for s = 1:numel(subjects)
    try
        subj = subjects(s);
        safeSubj = regexprep(subj,'[<>:"/\\|?*]+','_');

        % --- Load per-trial table (built earlier on the continuous stream)
        tsv = fullfile(inEpochs, safeSubj + "_trialtable.csv");
        if ~isfile(tsv)
            warning('[%s] trialtable missing → skipping subject.', safeSubj);
            continue
        end
        T = readtable(tsv, 'Delimiter',';', 'TextType','string');
        % Ensure the key columns exist
        if ~ismember("cue_urev", string(T.Properties.VariableNames)) || ...
           ~ismember("targ_urev", string(T.Properties.VariableNames))
            warning('[%s] trialtable lacks cue_urev/targ_urev → skipping subject.', safeSubj);
            continue
        end
        % Normalize error_type/resp_side/targ_side to string columns
        if ismember("error_type", string(T.Properties.VariableNames))
            T.error_type = string(T.error_type);
        else
            T.error_type = strings(height(T),1);
        end
        if ismember("resp_side", string(T.Properties.VariableNames))
            T.resp_side = upper(string(T.resp_side));
        else
            T.resp_side = strings(height(T),1);
        end
        if ismember("targ_side", string(T.Properties.VariableNames))
            T.targ_side = upper(string(T.targ_side));
        else
            T.targ_side = strings(height(T),1);
        end

        if ~ismember("rt_ms", string(T.Properties.VariableNames))
        T.rt_ms = nan(height(T),1);
        end


        %% ---------------- CUE (catch-only): drop FALSE ALARMS ----------------
        fcue = fullfile(inEpochs, safeSubj + "_cue_clean.set");
        if isfile(fcue)
        EEG = pop_loadset('filename', char(safeSubj + "_cue_clean.set"),  'filepath', char(inEpochs));           
        m = EEG.trials;
            keep = true(1, m);
            assert(numel(keep) == EEG.trials, 'keep vector length mismatch with EEG.trials');

            for k = 1:m
                % Get event at time 0 in this epoch
                evIdx = EEG.epoch(k).event;
                lat   = EEG.epoch(k).eventlatency;
                if iscell(lat),   lat = cellfun(@(x) x, lat); end
                if iscell(evIdx), evIdx = cellfun(@(x) x, evIdx); end
                z = find(abs(lat) <= epsZero, 1);  % cue at 0 ms
                if isempty(z), continue; end

                % Prefer epoch.eventurevent if present; else go via EEG.event
                u = NaN;
                if isfield(EEG.epoch, 'eventurevent') && ~isempty(EEG.epoch(k).eventurevent)
                    ue = EEG.epoch(k).eventurevent;
                    if iscell(ue), ue = cellfun(@(x) x, ue); end
                    if numel(ue) >= z && ~isempty(ue(z)) && all(isfinite(ue(z)))
                        u = double(ue(z));
                    end
                end
                if ~(isfinite(u))
                    e = evIdx(z);
                    if iscell(e), e = e{1}; end
                    if ~(isscalar(e) && e>=1 && e<=numel(EEG.event))
                        keep(k) = false;  % conservative: drop if the index is invalid
                        continue
                    end
                    if isfield(EEG.event,'urevent') && ~isempty(EEG.event(e).urevent)
                        u = double(EEG.event(e).urevent);
                    end
                end
                if ~(isfinite(u)), continue; end

                % Look up this cue in trial table and drop if false alarm
                ii = find(T.cue_urev == u, 1);
                if ~isempty(ii)
                    et = lower(strtrim(T.error_type(ii)));
                    if et == "false_alarm"
                        keep(k) = false;
                    end
                end
            end

            nDrop = sum(~keep);
            fprintf('[%s|CUE] Removing %d/%d false-alarm epochs.\n', safeSubj, nDrop, m);

            % Save audit CSV
            behCSV = fullfile(outEpochsBeh, sprintf('%s_cue_behavior_kept.csv', safeSubj));
            writetable(table((1:m).', keep.', 'VariableNames', {'epoch','kept'}), behCSV, 'Delimiter',';');

            % Apply and save
            EEG.reject.rejmanual  = ~keep;
            EEG.reject.rejmanualE = repmat(~keep, EEG.nbchan, 1);
            if any(~keep)
                EEG = pop_rejepoch(EEG, EEG.reject.rejmanual, 0);
            end
            EEG.setname = char(safeSubj + "_cue_clean_beh");
            pop_saveset(EEG, 'filename', char(EEG.setname + ".set"), 'filepath', char(outEpochsBeh));
            fprintf('  → saved %s\n', fullfile(outEpochsBeh, EEG.setname + ".set"));
            clear EEG  % <— optional memory helper
        else
            fprintf('[%s|CUE] No cue_clean set—skipped.\n', safeSubj);
        end

%% ---------------- TARGET: drop WRONG-SIDE and RT out-of-window (also drop misses) ----------------
ftarg = fullfile(inEpochs, safeSubj + "_targ_clean.set");
if isfile(ftarg)
    EEG = pop_loadset('filename', char(safeSubj + "_targ_clean.set"), 'filepath', char(inEpochs));
    m = EEG.trials;
    keep = true(1, m);
    assert(numel(keep) == EEG.trials, 'keep vector length mismatch with EEG.trials');

    % RT window (ms)
    rt_min_ms = 200;
    rt_max_ms = 1000;

    % Prepare an audit table that will be filled with the loop
    epoch_idx   = (1:m).';
    targ_urev_v = nan(m,1);

    % initialize to empty strings "" (NOT missing)
    errtype_v   = repmat("", m, 1);
    resp_side_v = repmat("", m, 1);
    targ_side_v = repmat("", m, 1);
    reason_v    = repmat("", m, 1);

    rt_ms_v     = nan(m,1);


    for k = 1:m
        % ----- find the target-at-zero event for this epoch -----
        evIdx = EEG.epoch(k).event;
        lat   = EEG.epoch(k).eventlatency;
        if iscell(lat),   lat = cellfun(@(x) x, lat); end
        if iscell(evIdx), evIdx = cellfun(@(x) x, evIdx); end
        z = find(abs(lat) <= epsZero, 1);
        if isempty(z)      % shouldn't happen; mark to drop conservatively
            keep(k)     = false;
            reason_v(k) = "no_zero_event";
            continue
        end

        % Robust urevent for this epoch's target
        u = NaN;
        if isfield(EEG.epoch, 'eventurevent') && ~isempty(EEG.epoch(k).eventurevent)
            ue = EEG.epoch(k).eventurevent;
            if iscell(ue), ue = cellfun(@(x) x, ue); end
            if numel(ue) >= z && ~isempty(ue(z)) && all(isfinite(ue(z)))
                u = double(ue(z));
            end
        end
        if ~(isfinite(u))
            e = evIdx(z); if iscell(e), e = e{1}; end
        if ~(isscalar(e) && e>=1 && e<=numel(EEG.event))
            keep(k)     = false;
            % choose a reason field that exists in this block:
            if exist('reason_v','var'); reason_v(k) = "bad_event_index"; end
            continue
        end

            if isfield(EEG.event,'urevent') && ~isempty(EEG.event(e).urevent)
                u = double(EEG.event(e).urevent);
            end
        end
        targ_urev_v(k) = u;

        % Match to trial table row
        ii = [];
        if isfinite(u)
            ii = find(T.targ_urev == u, 1);
        end
        if isempty(ii)
            keep(k)     = false;
            reason_v(k) = "no_trialtable_match";
            continue
        end

        % Pull behavior fields
        et  = lower(strtrim(T.error_type(ii)));    % "correct"|"wrong_side"|"anticipation"|"late"|"miss"|...
        rs0 = upper(strtrim(T.resp_side(ii)));     % "L"|"R"|"" 
        ts0 = upper(strtrim(T.targ_side(ii)));     % "L"|"R"|"" 
        rtm = T.rt_ms(ii);                         % numeric (NaN for miss)

        % Store into audit vectors (avoid <missing>)
        if ~ismissing(et),  errtype_v(k)   = string(et);  else, errtype_v(k)   = ""; end
        if ~ismissing(rs0), resp_side_v(k) = string(rs0); else, resp_side_v(k) = ""; end
        if ~ismissing(ts0), targ_side_v(k) = string(ts0); else, targ_side_v(k) = ""; end
        rt_ms_v(k) = rtm;

        % --------- decision logic ----------
        % (A) wrong side (explicit or inferred)
        wrong_side = (errtype_v(k) == "wrong_side");
        if ~wrong_side
            rs = resp_side_v(k); ts = targ_side_v(k);
            if (rs=="L" && ts=="R") || (rs=="R" && ts=="L")
                wrong_side = true;
            end
        end

        % (B) RT window
        has_rt  = isfinite(rtm) && ~isnan(rtm);
        too_fast = has_rt && (rtm < rt_min_ms);
        too_slow = has_rt && (rtm > rt_max_ms);

        % (C) miss (no RT)
        is_miss = (errtype_v(k) == "miss") || ~has_rt;

        % We KEEP only epochs that are correct AND in the RT window
        is_correct = (errtype_v(k) == "correct");
        in_window  = has_rt && ~too_fast && ~too_slow;

        keep(k) = is_correct && in_window;

        if ~keep(k)
            % Build reasons from tokens (avoid leading ';' and missing)
            reasons = {};
            if wrong_side, reasons{end+1} = 'wrong_side'; end %#ok<AGROW>
            if is_miss,    reasons{end+1} = 'miss/no_rt'; end %#ok<AGROW>
            if too_fast,   reasons{end+1} = sprintf('rt<%d', rt_min_ms); end %#ok<AGROW>
            if too_slow,   reasons{end+1} = sprintf('rt>%d', rt_max_ms); end %#ok<AGROW>
            if isempty(reasons)
                % Fallback to error_type (already lower-case)
                reasons = { char(errtype_v(k)) };
            end
            reason_v(k) = string(strjoin(reasons, ';'));
        else
            reason_v(k) = "";  % kept
        end
    end

    % ----- Console summary + list of dropped epochs -----
    nDrop = sum(~keep);
    fprintf('[%s|TARG] Removing %d/%d epochs (reasons: wrong-side and/or RT out-of-window and/or miss).\n', safeSubj, nDrop, m);

    if nDrop > 0
        dropped_idx = find(~keep);
        fprintf('  Dropped epochs (epoch  targ_urev  reason  rt_ms  resp→targ):\n');
        max_print = min(40, numel(dropped_idx));  % print up to 40 lines; see CSV for full list
        for q = 1:max_print
            e = dropped_idx(q);

            % target urevent as text
            tu = targ_urev_v(e);
            if isfinite(tu), tu_str = num2str(tu);
            else,            tu_str = 'NaN'; end

            % rt as text
            rt = rt_ms_v(e);
            if isfinite(rt), rt_str = sprintf('%.1f', rt);
            else,            rt_str = 'NaN'; end

            % safe strings (avoid <missing> -> '')
            rstr = reason_v(e); if ismissing(rstr), rstr = ""; end
            rs   = resp_side_v(e); if ismissing(rs), rs = ""; end
            ts   = targ_side_v(e); if ismissing(ts), ts = ""; end

            fprintf('   %4d   %7s   %s   %s   %s→%s\n', ...
                e, tu_str, char(rstr), rt_str, char(rs), char(ts));
        end
        if numel(dropped_idx) > max_print
            fprintf('   ... and %d more (see CSV)\n', numel(dropped_idx)-max_print);
        end
    end

    % ----- Write a detailed audit CSV -----
    behCSV = fullfile(outEpochsBeh, sprintf('%s_target_behavior_kept.csv', safeSubj));
    Tbeh = table( ...
        epoch_idx, keep.', targ_urev_v, errtype_v, resp_side_v, targ_side_v, rt_ms_v, reason_v, ...
        'VariableNames', {'epoch','kept','targ_urevent','error_type','resp_side','targ_side','rt_ms','drop_reason'});
    writetable(Tbeh, char(behCSV), 'Delimiter',';');

    % ----- Apply and save -----
    EEG.reject.rejmanual  = ~keep;
    EEG.reject.rejmanualE = repmat(~keep, EEG.nbchan, 1);
    if any(~keep)
        EEG = pop_rejepoch(EEG, EEG.reject.rejmanual, 0);
    end
    EEG.setname = char(safeSubj + "_targ_clean_beh");
    pop_saveset(EEG, 'filename', char(EEG.setname + ".set"), 'filepath', char(outEpochsBeh));
    fprintf('  → saved %s\n', fullfile(outEpochsBeh, EEG.setname + ".set"));
    clear EEG  % <— optional memory helper
else
    fprintf('[%s|TARG] No targ_clean set—skipped.\n', safeSubj);
end



    catch ME
        warning('[BehaviorFilter] %s failed: %s - %s', safeSubj, ME.identifier, ME.message);
        continue
    end
end

fprintf('[BehaviorFilter] Done. Files written to: %s\n', outEpochsBeh);





%% ================== Cohort summary: epoched vs. kept vs. dropped ==================
% Folders
inEpochsBeh = '...\BehaviorClean';
[inEpochs, ~, ~] = fileparts(inEpochsBeh);  % parent folder that holds *_clean.set

% EEGLAB init (adjust path if needed)
addpath('...\eeglab2025.0.0');
eeglab; close;  % init quietly

% Find post-behavior sets and build subject list
filesCueB  = dir(fullfile(inEpochsBeh, '*_cue_clean_beh.set'));
filesTargB = dir(fullfile(inEpochsBeh, '*_targ_clean_beh.set'));

subs = strings(0,1);
for k = 1:numel(filesCueB)
    nm = erase(string(filesCueB(k).name), '_cue_clean_beh.set');
    subs(end+1) = nm; %#ok<SAGROW>
end
for k = 1:numel(filesTargB)
    nm = erase(string(filesTargB(k).name), '_targ_clean_beh.set');
    subs(end+1) = nm; %#ok<SAGROW>
end
subs = unique(subs(subs~=""));

% Storage
Subject        = strings(0,1);
nCueBefore     = []; nCueAfter  = []; nCueDrop  = [];
nTargBefore    = []; nTargAfter = []; nTargDrop = [];

fprintf('\n[Behavior cohort summary] Epoched (before) vs Kept (after) vs Dropped\n');
fprintf('%-30s | %7s %7s %7s || %7s %7s %7s\n', ...
    'Subject', 'Cue_B', 'Cue_A', 'Cue_D', 'Targ_B', 'Targ_A', 'Targ_D');
fprintf(repmat('-',1,30)); fprintf('+'); fprintf(repmat('-',1,26)); fprintf('+'); fprintf(repmat('-',1,26)); fprintf('\n');

for i = 1:numel(subs)
    subj = subs(i);

    % ---------- Cue ----------
    nB_cue = NaN; nA_cue = NaN;
    fCueBefore = fullfile(inEpochs,    subj + "_cue_clean.set");
    fCueAfter  = fullfile(inEpochsBeh, subj + "_cue_clean_beh.set");

    if isfile(fCueBefore)
        EEGcB = pop_loadset('filename', char(subj + "_cue_clean.set"), ...
                            'filepath', char(inEpochs), 'loadmode', 'info');
        nB_cue = EEGcB.trials;
    end
    if isfile(fCueAfter)
        EEGcA = pop_loadset('filename', char(subj + "_cue_clean_beh.set"), ...
                            'filepath', char(inEpochsBeh), 'loadmode', 'info');
        nA_cue = EEGcA.trials;
    end
    nD_cue = nB_cue - nA_cue;  % will be NaN if any side missing

    % ---------- Target ----------
    nB_targ = NaN; nA_targ = NaN;
    fTargBefore = fullfile(inEpochs,    subj + "_targ_clean.set");
    fTargAfter  = fullfile(inEpochsBeh, subj + "_targ_clean_beh.set");

    if isfile(fTargBefore)
        EEGtB = pop_loadset('filename', char(subj + "_targ_clean.set"), ...
                            'filepath', char(inEpochs), 'loadmode', 'info');
        nB_targ = EEGtB.trials;
    end
    if isfile(fTargAfter)
        EEGtA = pop_loadset('filename', char(subj + "_targ_clean_beh.set"), ...
                            'filepath', char(inEpochsBeh), 'loadmode', 'info');
        nA_targ = EEGtA.trials;
    end
    nD_targ = nB_targ - nA_targ;

    % Print row
    fprintf('%-30s | %7s %7s %7s || %7s %7s %7s\n', ...
        subj, mat2str(nB_cue), mat2str(nA_cue), mat2str(nD_cue), ...
              mat2str(nB_targ), mat2str(nA_targ), mat2str(nD_targ));

    % Collect row
    Subject(end+1,1)     = subj; %#ok<SAGROW>
    nCueBefore(end+1,1)  = nB_cue; %#ok<SAGROW>
    nCueAfter(end+1,1)   = nA_cue; %#ok<SAGROW>
    nCueDrop(end+1,1)    = nD_cue; %#ok<SAGROW>
    nTargBefore(end+1,1) = nB_targ; %#ok<SAGROW>
    nTargAfter(end+1,1)  = nA_targ; %#ok<SAGROW>
    nTargDrop(end+1,1)   = nD_targ; %#ok<SAGROW>
end

fprintf(repmat('-',1,84)); fprintf('\nDone.\n');

% ------- Write single cohort CSV -------
Summary = table(Subject, ...
    nCueBefore,  nCueAfter,  nCueDrop, ...
    nTargBefore, nTargAfter, nTargDrop, ...
    'VariableNames', {'subject', ...
        'n_cue_epoched','n_cue_kept','n_cue_dropped', ...
        'n_targ_epoched','n_targ_kept','n_targ_dropped'});

outCSV = fullfile(inEpochsBeh, 'behavior_cohort_summary.csv');
writetable(Summary, outCSV, 'Delimiter',';');

fprintf('Wrote cohort summary: %s\n', outCSV);