%-----------------------------------------------------------------------%
% EEG time-frequency decomposition for the Posner spatial attention task
%-----------------------------------------------------------------------%
% uses: MATLAB R2025a | EEGLAB 2025.0.0 | FieldTrip

% Author: Leonie Polzer

% Description:
% - Batch time-frequency analysis of cleaned, target-locked EEG data
% - Uses FieldTrip Morlet wavelets (8–60 Hz) to extract oscillatory power,
%   inter-trial phase coherence (ITPC), and inter-site phase coherence (ISPC).
% - Calculates alpha, beta, and gamma measures across occipital, parietal,
%   and frontal ROIs, separately for ipsi- and contralateral hemispheres,
%   task conditions, target sides, and cue-target intervals (SOA).
% - Links behavioral data (reaction times and accuracy) to EEG trials.
% - Exports trial-wise and participant-level CSVs, time-frequency matrices,
%   and connectivity data for subsequent statistical analyses.


%% ---- Environment -----
% clear space
clearvars; close all; clc;

% eeglab
addpath('...\eeglab2025.0.0');
eeglab; set(0,'DefaultFigureWindowStyle','normal');

% FieldTrip 
addpath('...\fieldtrip-20250915');   
ft_defaults;

% Paths
inEpochs = '...\BehaviorClean';
outBase  = '...\TF_Final_fieldtrip';

outTrial = fullfile(outBase,'trialwise_fieldtrip');
outITC   = fullfile(outBase,'itc_subject_fieldtrip');
outISPC = fullfile(outBase,'ispc_subject_fieldtrip');
outDense = fullfile(outBase,'roi_dense_fieldtrip');
outConn  = fullfile(outBase,'conn_cache_fieldtrip');

if ~isfolder(outBase),  mkdir(outBase);  end
if ~isfolder(outTrial), mkdir(outTrial); end
if ~isfolder(outITC),   mkdir(outITC);   end
if ~isfolder(outISPC), mkdir(outISPC); end
if ~isfolder(outDense), mkdir(outDense); end
if ~isfolder(outConn),  mkdir(outConn);  end

diary(fullfile(outBase, sprintf('TF_fieldtrip_runlog_%s.txt', datestr(now,'yyyymmdd_HHMMSS'))));
diary on;

% File discovery
filesTarg = dir(fullfile(inEpochs, '*_targ_clean_beh.set'));
fprintf('[TF-FT] Found %d target-locked datasets.\n', numel(filesTarg));

%% Define settings and variables
%%% ---- Frequency grid and cycles ----
N_FREQ = 30;
freqs  = logspace(log10(8), log10(60), N_FREQ);

c_lo = 4;
c_hi = 9;
cycles = c_lo + (c_hi - c_lo) * ( log(freqs/8) / log(60/8) );
cycles = max(min(cycles, 9), 4);

n_f = numel(freqs);

% Snap nominal 5-ms targets to the true sampled axis; effective step is ~4/6 ms at 500 Hz
crop_tmin_ms  = -500;
crop_tmax_ms  = 1000;
crop_step_ms  = 5;

crop_times_ms_nominal = (crop_tmin_ms:crop_step_ms:crop_tmax_ms); 

% Baseline & condition field names 
BASE_T1      = 'bl_precue_t1';
BASE_T2      = 'bl_precue_t2';
SIDE_TARG_F  = 'targ_side';
COND_F       = 'validity';

%%% ---- Analysis bins ----
fine_ids = 1:4;
fine_win = [ ...
    0   125; ...   % [0,125)
    125 250; ...   % [125,250)
    250 375; ...   % [250,375)
    375 501  ...   % [375,501) effectively <=500
];

coll_ids = [101 102 999];
coll_win = [ ...
    0   250; ...   % [0,250)
    250 501; ...   % [250,501)
    0   501  ...   % [0,501)  -> FULL 0–500 ms average
];


%%% ---- Bands (Hz) ----
band_names  = {'alpha','beta','gamma_bb'};
band_bounds = [ 8 12; 14 29; 30 60 ];
notch_excl  = [48 52];

band_sel = struct();
for bb = 1:numel(band_names)
    f_lo_b = band_bounds(bb,1); f_hi_b = band_bounds(bb,2);
    keep = (freqs >= f_lo_b) & (freqs <= f_hi_b);
    keep = keep & ~((freqs >= notch_excl(1)) & (freqs <= notch_excl(2)));
    band_sel.(band_names{bb}) = find(keep);
end
n_freq_in_band = struct('alpha',numel(band_sel.alpha), 'beta',numel(band_sel.beta), 'gamma_bb',numel(band_sel.gamma_bb));

%%% ---- ROIs ----
ROI.occi.all = {'PO7','PO8','PO3','PO4','O1','O2'};
ROI.occi.L   = {'PO7','PO3','O1'};
ROI.occi.R   = {'PO8','PO4','O2'};

ROI.pari.all = {'P5','P6','P1','P2','P3','P4'};
ROI.pari.L   = {'P5','P1','P3'};
ROI.pari.R   = {'P6','P2','P4'};

ROI.front.all = {'FC3','FC4','FC1','FC2','F1','F2'};
ROI.front.L   = {'FC3','FC1','F1'};
ROI.front.R   = {'FC4','FC2','F2'};

roi_names = {'occi','pari','front'};

% Safety thresholds / assertions 
min_pre_ms  = 1000;
min_post_ms = 1000;

% CSV headers
TW_hdr = {'subj_id','setname','epoch_idx','urevent','which','validity','side','SOA_ms','SOA_label', ...
          'ROI','hemi','band','time_window','tbin_start_ms','tbin_end_ms', ...
          'dB_total','RT_ms','acc','n_elec_hemi','n_freq_in_band','n_time_in_bin'};

ITC_hdr = {'subj_id','which','validity','side','SOA_label','ROI','hemi','band','time_window', ...
           'tbin_start_ms','tbin_end_ms','ITC','ITPCz','n_trials_used','n_freq_in_band','n_time_in_bin'};

ISPC_hdr = {'subj_id','which','roi1','roi2','laterality','validity','side','SOA_label', ...
            'band','time_window','tbin_start_ms','tbin_end_ms', ...
            'ISPC','ISPCz','n_trials','n_pairs','n_freq_in_band','n_time_in_bin'};

summary_rows = {};
summary_hdr  = {'subject','n_trl_targ','n_cued_targ','n_uncued_targ','behavior_mapped','rt_nonNaN','acc_nonNaN'};

fprintf('[TF-FT] Setup complete. Starting TARGET-only batch…\n');

%%
 % % ---------------------------
 % % SUBJECT LOOP (TARGET ONLY)
 % % ---------------------------

for iFile = 1:numel(filesTarg)

    targ_name   = filesTarg(iFile).name;
    subj_prefix = regexprep(targ_name, '_targ_clean_beh\.set$', '');
    subj_id     = subj_prefix;

    fprintf('\n[%s] [SUBJ %d/%d] Starting %s\n', datestr(now,'HH:MM:SS'), iFile, numel(filesTarg), targ_name);
    t_subj = tic;

    try
        % ---------- Load target-locked dataset ----------
        EEG = pop_loadset('filename', targ_name, 'filepath', inEpochs);
        EEG = eeg_checkset(EEG, 'eventconsistency');

        whichstr = 'targ';
        fs    = EEG.srate;
        nChan = EEG.nbchan;
        nTr   = EEG.trials;
        tvec  = EEG.times;

        if (abs(min(tvec)) < min_pre_ms) || (max(tvec) < min_post_ms)
            error('Target epochs too short: need ≥%d ms pre and ≥%d ms post. Found [%g %g] ms.', ...
                  min_pre_ms, min_post_ms, min(tvec), max(tvec));
        end

        % ---------- ROI channel indices 
        all_labels       = string({EEG.chanlocs.labels});
        all_labels_lower = lower(all_labels);

        roi_idx = struct();
        warn_missing = {};
        for r = 1:numel(roi_names)
            rn = roi_names{r};
            roi_idx.(rn).all = [];
            roi_idx.(rn).L   = [];
            roi_idx.(rn).R   = [];
            sets = {'all','L','R'};
            for s = 1:numel(sets)
                list = ROI.(rn).(sets{s});
                idx_here = [];
                for k = 1:numel(list)
                    hit = find(all_labels_lower == lower(string(list{k})), 1, 'first');
                    if ~isempty(hit)
                        idx_here(end+1) = hit; %#ok<AGROW>
                    else
                        warn_missing{end+1} = sprintf('%s:%s channel "%s"', rn, sets{s}, list{k}); %#ok<AGROW>
                    end
                end
                roi_idx.(rn).(sets{s}) = unique(idx_here);
            end
        end
        if ~isempty(warn_missing)
            fprintf(2, '[TF-FT][WARN] Missing/misnamed ROI channels (target):\n');
            for w = 1:numel(warn_missing), fprintf(2, '  - %s\n', warn_missing{w}); end
        end

        % Nominal per-hemi channel counts
        n_elec_hemi_nominal = struct();
        for r = 1:numel(roi_names)
            rn = roi_names{r};
            n_elec_hemi_nominal.(rn) = struct();
            n_elec_hemi_nominal.(rn).ipsi   = struct('L', numel(roi_idx.(rn).L), 'R', numel(roi_idx.(rn).R));
            n_elec_hemi_nominal.(rn).contra = struct('L', numel(roi_idx.(rn).R), 'R', numel(roi_idx.(rn).L));
        end

        % ---------- Per-trial metadata (target) ----------
        side_raw_targ = repmat("R", 1, nTr);
        if isfield(EEG,'epoch') && isfield(EEG.epoch,SIDE_TARG_F)
            try
                side_raw_targ = string({EEG.epoch.(SIDE_TARG_F)});
                if numel(side_raw_targ) ~= nTr
                    tmp = repmat("R", 1, nTr);
                    nCopy = min(numel(side_raw_targ), nTr);
                    tmp(1:nCopy) = side_raw_targ(1:nCopy);
                    side_raw_targ = tmp;
                end
            catch
                side_raw_targ = repmat("R", 1, nTr);
            end
        elseif isfield(EEG,'epoch') && isfield(EEG.epoch,'side')
            try
                side_raw_targ = string({EEG.epoch.side});
            catch
                side_raw_targ = repmat("R", 1, nTr);
            end
        end
        side_raw_targ = upper(strtrim(side_raw_targ));

        % ---------- SIDE sanity check ----------
        side_norm = upper(strtrim(side_raw_targ));
        ok_side = (side_norm=="L") | (side_norm=="R");
        n_bad_side = sum(~ok_side);

        if n_bad_side > 0
            fprintf(2,'[TF-FT][WARN][%s] %d/%d trials have invalid targ_side labels (not L/R). Examples:\n', ...
                subj_id, n_bad_side, nTr);
            bad_idx = find(~ok_side);
            for ii = 1:min(10,numel(bad_idx))
                trb = bad_idx(ii);
                fprintf(2,'  trial %d: "%s"\n', trb, char(side_raw_targ(trb)));
            end
        end

        validity_raw_targ = repmat("", 1, nTr);
        if isfield(EEG,'epoch') && isfield(EEG.epoch,COND_F)
            try
                validity_raw_targ = string({EEG.epoch.(COND_F)});
                if numel(validity_raw_targ) ~= nTr
                    tmp = repmat("", 1, nTr);
                    nCopy = min(numel(validity_raw_targ), nTr);
                    tmp(1:nCopy) = validity_raw_targ(1:nCopy);
                    validity_raw_targ = tmp;
                end
            catch
                validity_raw_targ = repmat("", 1, nTr);
            end
        end

        % SOA label: normalize empties to 'no_cue'
        SOA_label_targ = repmat("no_cue", 1, nTr);
        if isfield(EEG,'epoch') && isfield(EEG.epoch,'targ_soa')
            try
                tmpSOA = string({EEG.epoch.targ_soa});
                if numel(tmpSOA) ~= nTr
                    tmp = repmat("no_cue", 1, nTr);
                    nCopy = min(numel(tmpSOA), nTr);
                    tmp(1:nCopy) = tmpSOA(1:nCopy);
                    tmpSOA = tmp;
                end
                for k = 1:nTr
                    if strlength(strtrim(tmpSOA(k)))==0
                        SOA_label_targ(k) = "no_cue";
                    else
                        SOA_label_targ(k) = lower(strtrim(tmpSOA(k)));
                    end
                end
            catch
                SOA_label_targ = repmat("no_cue", 1, nTr);
            end
        end

        SOA_raw_targ = nan(1, nTr);
        if isfield(EEG,'epoch') && isfield(EEG.epoch,'targ_soa_ms')
            try
                SOA_raw_targ = [EEG.epoch.targ_soa_ms];
                if numel(SOA_raw_targ) ~= nTr
                    tmp = nan(1, nTr); nCopy = min(numel(SOA_raw_targ), nTr);
                    tmp(1:nCopy) = SOA_raw_targ(1:nCopy); SOA_raw_targ = tmp;
                end
            catch
                SOA_raw_targ = nan(1, nTr);
            end
        else
            for tr = 1:nTr
                if SOA_label_targ(tr)=="short", SOA_raw_targ(tr)=100;
                elseif SOA_label_targ(tr)=="long", SOA_raw_targ(tr)=250;
                else, SOA_raw_targ(tr)=NaN;
                end
            end
        end

        % validity 'no_cue' -> SOA_label 'no_cue'
        is_nc = strcmpi(validity_raw_targ,'no_cue');
        SOA_label_targ(is_nc) = "no_cue";

        

        % ---------- Behavior merge ----------
        RT_ms_targ   = nan(1, nTr);
        acc_targ     = nan(1, nTr);
        urevent_targ = nan(1, nTr);

        Tbeh = [];
        if isfield(EEG,'etc') && isfield(EEG.etc,'trialtable') && ~isempty(EEG.etc.trialtable)
            tt = EEG.etc.trialtable;
            if isstruct(tt) && isfield(tt,'trialtable') && ~istable(tt)
                tt = tt.trialtable;
            end
            if istable(tt)
                Tbeh = tt;
            elseif istimetable(tt)
                Tbeh = timetable2table(tt,'ConvertRowTimes',true);
            elseif isstruct(tt)
                try
                    Tbeh = struct2table(tt);
                catch
                    Tbeh = [];
                end
            elseif iscell(tt) && numel(tt)==1 && istable(tt{1})
                Tbeh = tt{1};
            end
        end

        if isempty(Tbeh) || ~istable(Tbeh) || height(Tbeh)==0
            fprintf(2,'[TF-FT][BEH][%s] trialtable missing/unusable. RT/ACC remain NaN.\n', subj_id);
        else
            behVarsOrig  = string(Tbeh.Properties.VariableNames);
            behVarsLower = lower(behVarsOrig);

            key_idx = [];
            for cand = ["targ_urev","targ_urevent","target_urev"]
                idx = find(behVarsLower==cand, 1);
                if ~isempty(idx), key_idx = idx; break; end
            end

            if isempty(key_idx)
                fprintf(2,'[TF-FT][BEH][%s] trialtable lacks target key column.\n', subj_id);
            else
                key_raw = Tbeh.(behVarsOrig(key_idx));

                if isnumeric(key_raw)
                    key_vec_all = double(key_raw);
                elseif islogical(key_raw)
                    key_vec_all = double(key_raw);
                else
                    try
                        s = string(key_raw);
                    catch
                        s = strings(size(key_raw));
                        for ii=1:numel(key_raw)
                            try, s(ii) = string(key_raw(ii)); catch, s(ii) = ""; end
                        end
                    end
                    s = strrep(s, ',', '.');
                    key_vec_all = str2double(s);
                end

                keep_rows = ~isnan(key_vec_all);
                key_vec   = key_vec_all(keep_rows);

                % RT
                has_rt = any(behVarsLower=="rt_ms") || any(behVarsLower=="rt") || ...
                         any(behVarsLower=="reaction_time_ms") || any(behVarsLower=="reaction_time");
                if has_rt
                    if     any(behVarsLower=="rt_ms"),             rt_idx = find(behVarsLower=="rt_ms",1);
                    elseif any(behVarsLower=="reaction_time_ms"),  rt_idx = find(behVarsLower=="reaction_time_ms",1);
                    elseif any(behVarsLower=="reaction_time"),     rt_idx = find(behVarsLower=="reaction_time",1);
                    else                                           rt_idx = find(behVarsLower=="rt",1);
                    end
                    vr = Tbeh.(behVarsOrig(rt_idx));
                    if isnumeric(vr)
                        val_rt = double(vr);
                    elseif islogical(vr)
                        val_rt = double(vr);
                    else
                        try, s = string(vr); catch, s = strings(size(vr)); end
                        s = strrep(s,',','.');
                        val_rt = str2double(s);
                    end
                    val_rt = val_rt(keep_rows);
                else
                    val_rt = [];
                end

                % ACC
                has_acc = any(behVarsLower=="acc") || any(behVarsLower=="accuracy") || ...
                          any(behVarsLower=="correct") || any(behVarsLower=="is_correct");
                if has_acc
                    if     any(behVarsLower=="acc"),       acc_idx = find(behVarsLower=="acc",1);
                    elseif any(behVarsLower=="accuracy"),  acc_idx = find(behVarsLower=="accuracy",1);
                    elseif any(behVarsLower=="correct"),   acc_idx = find(behVarsLower=="correct",1);
                    else                                   acc_idx = find(behVarsLower=="is_correct",1);
                    end
                    va = Tbeh.(behVarsOrig(acc_idx));
                    if islogical(va)
                        val_acc = double(va);
                    elseif isnumeric(va)
                        val_acc = double(va ~= 0);
                    else
                        try, s = string(va); catch, s = strings(size(va)); end
                        s_low = lower(strtrim(s));
                        val_acc = nan(size(s_low));
                        pos = ismember(s_low, ["1","true","wahr","yes","ja","correct","korrekt","richtig","ok"]);
                        neg = ismember(s_low, ["0","false","falsch","no","nein","incorrect","inkorrekt","nicht korrekt","wrong"]);
                        val_acc(pos) = 1;  val_acc(neg) = 0;
                        m = isnan(val_acc);
                        if any(m)
                            s_num = strrep(s_low(m),',','.');
                            val_acc(m) = double(str2double(s_num) ~= 0);
                        end
                    end
                    val_acc = val_acc(keep_rows);
                else
                    val_acc = [];
                end

                if isempty(key_vec)
                    fprintf(2,'[TF-FT][BEH][%s] All target keys NaN — no behavior mapped.\n', subj_id);
                else
                    Bset = unique(key_vec(:).');
                    if ~isempty(val_rt)
                        Mrt = containers.Map(num2cell(key_vec(:).'), num2cell(val_rt(:).'));
                        has_rt = true;
                    else
                        Mrt = containers.Map('KeyType','double','ValueType','double');
                        has_rt = false;
                    end
                    if ~isempty(val_acc)
                        Macc = containers.Map(num2cell(key_vec(:).'), num2cell(val_acc(:).'));
                    else
                        Macc = containers.Map('KeyType','double','ValueType','double');
                    end

                    epsZero = 1e-3;
                    for tr = 1:nTr
                        urevent_targ(tr) = NaN;
                        RT_ms_targ(tr)   = NaN;
                        acc_targ(tr)     = NaN;

                        if ~isfield(EEG,'epoch') || numel(EEG.epoch) < tr, continue; end
                        if ~isfield(EEG.epoch(tr),'event') || ~isfield(EEG.epoch(tr),'eventlatency'), continue; end

                        evIdx = EEG.epoch(tr).event;
                        lat   = EEG.epoch(tr).eventlatency;
                        if iscell(lat),   lat   = cellfun(@double, lat);   end
                        if iscell(evIdx), evIdx = cellfun(@double, evIdx); end
                        if isempty(evIdx) || isempty(lat), continue; end

                        chosen_ev = NaN;
                        z = find(abs(lat) <= epsZero, 1, 'first');
                        if ~isempty(z)
                            e = evIdx(z);
                            if isfinite(e) && e>=1 && e<=numel(EEG.event)
                                chosen_ev = e;
                            end
                        end
                        if ~isfinite(chosen_ev)
                            [~,kmin] = min(abs(lat));
                            e = evIdx(kmin);
                            if isfinite(e) && e>=1 && e<=numel(EEG.event)
                                chosen_ev = e;
                            end
                        end

                        if isfinite(chosen_ev)
                            chosen_urev = double(EEG.event(chosen_ev).urevent);
                            if isfinite(chosen_urev) && any(chosen_urev == Bset)
                                urevent_targ(tr) = chosen_urev;
                                if has_rt && isKey(Mrt, chosen_urev),  RT_ms_targ(tr) = Mrt(chosen_urev); end
                                if         isKey(Macc, chosen_urev), acc_targ(tr)   = Macc(chosen_urev); end
                            end
                        end
                    end

                    n_map = sum(~isnan(urevent_targ));
                    fprintf('[VAL][%s] Behavior merge mapped %d/%d target epochs (%.1f%%)\n', ...
                        subj_id, n_map, nTr, 100*n_map/max(1,nTr));

                    coverage = n_map / max(1,nTr);

                    if coverage < 0.90
                        fprintf(2,'[TF-FT][BEH][WARNING][%s] Only %.1f%% of trials matched to behavior (%d/%d).\n', ...
                            subj_id, 100*coverage, n_map, nTr);
                    end
                end
            end
        end

        isLeft_targ = false(1, nTr);
        for tr = 1:nTr
            isLeft_targ(tr) = (side_raw_targ(tr)=="L");
        end

        % ---------- SIDE mapping consistency check ----------
        mis_L = find(strcmpi(side_raw_targ,"L") & ~isLeft_targ);
        mis_R = find(strcmpi(side_raw_targ,"R") &  isLeft_targ);
        if ~isempty(mis_L) || ~isempty(mis_R)
            fprintf(2,'[TF-FT][WARN][%s] isLeft_targ mismatch with side_raw_targ. mis_L=%d, mis_R=%d\n', ...
                subj_id, numel(mis_L), numel(mis_R));
        end

        fprintf('[CHK][%s] RT non-NaN: %d/%d | ACC non-NaN: %d/%d\n', ...
            subj_id, sum(~isnan(RT_ms_targ)), nTr, sum(~isnan(acc_targ)), nTr);

        % =================================================================
        % FieldTrip TF (Morlet wavelets) + EXPLICIT time mapping to crop grid
        % =================================================================
        ft_data = eeglab2fieldtrip(EEG, 'preprocessing', 'none');


       % --- Ensure FT channel order matches EEGLAB channel order ---
        ft_labels = string(ft_data.label(:));
        eeg_labels = string({EEG.chanlocs.labels});
        ft_lower = lower(strtrim(ft_labels));
        eeg_lower = lower(strtrim(eeg_labels(:)));
        
        if numel(ft_lower) ~= numel(eeg_lower)
            error('[TF-FT] Channel count mismatch EEGLAB vs FieldTrip (%d vs %d).', numel(eeg_lower), numel(ft_lower));
        end
        
        map_ft_to_eeg = zeros(numel(eeg_lower),1);
        for cc = 1:numel(eeg_lower)
            hit = find(ft_lower == eeg_lower(cc), 1, 'first');
            if isempty(hit)
                error('[TF-FT] Channel "%s" not found in FieldTrip labels.', eeg_labels{cc});
            end
            map_ft_to_eeg(cc) = hit;
        end
        

        % Ensure sampleinfo exists
        if ~isfield(ft_data, 'sampleinfo') || isempty(ft_data.sampleinfo)
            nTr_ft = numel(ft_data.trial);
            nS     = numel(ft_data.time{1});
            ft_data.sampleinfo = zeros(nTr_ft, 2);
            for tr = 1:nTr_ft
                a = (tr-1)*nS + 1;
                b = tr*nS;
                ft_data.sampleinfo(tr,:) = [a b];
            end
        end

        % Padding
        if ~isempty(ft_data.time) && ~isempty(ft_data.time{1})
            trial_len_s = double(ft_data.time{1}(end) - ft_data.time{1}(1));
        else
            trial_len_s = double((max(tvec)-min(tvec))/1000);
        end
        pad_s = max(ceil(trial_len_s), 2);

        cfg_tf = [];
        cfg_tf.method     = 'wavelet';
        cfg_tf.output     = 'fourier';
        cfg_tf.foi        = freqs;
        cfg_tf.width      = cycles;
        cfg_tf.toi        = ft_data.time{1};  % keep original sampling
        cfg_tf.pad        = pad_s;
        cfg_tf.channel    = 'all';
        cfg_tf.keeptrials = 'yes';

        fprintf('[%s][%s] Running FieldTrip TF: nTr=%d, nChan=%d, nFreq=%d, toi=%d samples, pad=%gs\n', ...
            datestr(now,'HH:MM:SS'), subj_id, EEG.trials, EEG.nbchan, numel(freqs), numel(cfg_tf.toi), cfg_tf.pad);

        t_tf = tic;
        freqFT = ft_freqanalysis(cfg_tf, ft_data);
        fprintf('[%s][%s] FieldTrip TF done in %.1fs\n', datestr(now,'HH:MM:SS'), subj_id, toc(t_tf));

        if ~isfield(freqFT,'fourierspctrm') || isempty(freqFT.fourierspctrm)
            error('[TF-FT] Missing fourierspctrm in FieldTrip output.');
        end

        % Robust dimensionality handling
        F = freqFT.fourierspctrm;
        szF = size(F);
        fprintf('[DBG][%s] fourierspctrm size = %s\n', subj_id, mat2str(szF));

        if numel(szF) == 5
            fprintf('[DBG][%s] fourierspctrm is 5-D -> averaging over dim5 (tapers)\n', subj_id);
            F = mean(F, 5, 'omitnan');
            szF = size(F);
            fprintf('[DBG][%s] fourierspctrm size after taper-mean = %s\n', subj_id, mat2str(szF));
        end

        if ndims(F) ~= 4
            error('[TF-FT] Unexpected fourierspctrm dimensionality: %d dims (size=%s).', ndims(F), mat2str(size(F)));
        end

        % Convert to complex single safely
        F_single = complex(single(real(F)), single(imag(F)));

        % Convert to cache-like order: [freq x chan x time x trial]
        Zfreq_all = permute(F_single, [3 2 4 1]);
        % --- Reorder channel dimension to EEGLAB order (NOW Zfreq_all exists) ---
        Zfreq_all = Zfreq_all(:, map_ft_to_eeg, :, :);

        % FieldTrip time/freq vectors
        times_ms_full = double(freqFT.time(:)' * 1000);
        freqs_ft      = double(freqFT.freq(:)');
        if any(diff(times_ms_full) <= 0)
            error('[TF-FT] freqFT.time is not strictly increasing. Cannot snap safely.');
        end

        fprintf('[DBG][%s] Zfreq_all size = %s\n', subj_id, mat2str(size(Zfreq_all)));
        fprintf('[DBG][%s] times_ms_full n=%d | freqs_ft n=%d\n', subj_id, numel(times_ms_full), numel(freqs_ft));

        % ---------- Frequency grid sanity check ----------
        if numel(freqs_ft) ~= numel(freqs)
            error('[TF-FT] FieldTrip returned %d freqs but requested %d.', numel(freqs_ft), numel(freqs));
        end

        max_abs_diff = max(abs(freqs_ft(:) - freqs(:)));
        if max_abs_diff > 1e-3
            fprintf(2,'[TF-FT][WARN][%s] freqs_ft differs from requested freqs (max |diff| = %.6f Hz). Using requested freqs for banding/metadata.\n', ...
                subj_id, max_abs_diff);
        end

        [md,ix] = max(abs(freqs_ft(:) - freqs(:)));
        if md > 1e-3
            fprintf(2,'[TF-FT][WARN][%s] freqs_ft differs (max |diff|=%.6f Hz) at idx=%d: ft=%.6f, req=%.6f\n', ...
                subj_id, md, ix, freqs_ft(ix), freqs(ix));
        end
        if md > 0.5
            error('[TF-FT] Frequency grid mismatch too large (%.3f Hz).', md);
        end

        % Snap nominal 5-ms targets to the TRUE sampled time axis,
        % then unique(...,'stable') so the effective step becomes 4/6 ms at 500 Hz.
        idx_crop_ft = zeros(1, numel(crop_times_ms_nominal));
        for ti = 1:numel(crop_times_ms_nominal)
            [~, ix] = min(abs(times_ms_full - crop_times_ms_nominal(ti)));
            idx_crop_ft(ti) = ix;
        end
        idx_crop_ft = unique(idx_crop_ft, 'stable');   % <-- THIS matches your old script
        
        crop_times_ms = times_ms_full(idx_crop_ft);    % <-- actual snapped times (often 4/6 ms steps)
        nTcrop        = numel(idx_crop_ft);
        
        % Crop to [freq x chan x nTcrop x trial]
        Zfreq_crop = Zfreq_all(:, :, idx_crop_ft, :);
        
        % From here on, "times_ms" is the ACTUAL snapped-to-sample grid
        times_ms = double(crop_times_ms);
        
        fprintf('[DBG][%s] crop_times_ms unique(diff) = %s | nTcrop=%d\n', subj_id, mat2str(unique(diff(times_ms))), nTcrop);
        
        % ---------- Build bins using HALF-OPEN intervals ----------
        % [lo, hi) except when hi == 501 (final inclusive window)
        
        fine_bins = cell(numel(fine_ids), 1);
        for b = 1:numel(fine_ids)
            lo = fine_win(b,1);
            hi = fine_win(b,2);
        
            if hi == 501
                fine_bins{b} = find(crop_times_ms >= lo & crop_times_ms <= 500);
            else
                fine_bins{b} = find(crop_times_ms >= lo & crop_times_ms < hi);
            end
        end
        
        coll_bins = cell(numel(coll_ids), 1);
        for b = 1:numel(coll_ids)
            lo = coll_win(b,1);
            hi = coll_win(b,2);
        
            if hi == 501
                coll_bins{b} = find(crop_times_ms >= lo & crop_times_ms <= 500);
            else
                coll_bins{b} = find(crop_times_ms >= lo & crop_times_ms < hi);
            end
        end

        fprintf('[DBG][%s] Zfreq_crop size = %s (expected approx [30 65 ~301 nTr])\n', subj_id, mat2str(size(Zfreq_crop)));

        % Sanity check: ensure 999 exists and covers 0..500
        ix999 = find(coll_ids == 999, 1, 'first');
        if isempty(ix999)
            error('[TF-FT] coll_ids missing 999 window.');
        end
        fprintf('[DBG][%s] n_time_in_bin (999, 0-500) = %d samples | t=[%.1f..%.1f] ms\n', ...
            subj_id, numel(coll_bins{ix999}), min(crop_times_ms(coll_bins{ix999})), max(crop_times_ms(coll_bins{ix999})));

        % ---------- Per-trial baseline indices on crop grid (per trial!) ----------
        if ~(isfield(EEG,'epoch') && isfield(EEG.epoch,BASE_T1) && isfield(EEG.epoch,BASE_T2))
            error('[TF-FT][BASE] Target epochs lack %s/%s.', BASE_T1, BASE_T2);
        end
        BL_targ_t1 = [EEG.epoch.(BASE_T1)];
        BL_targ_t2 = [EEG.epoch.(BASE_T2)];

        BL_idx_targ = cell(1, nTr);
        for tr = 1:nTr
            if isnan(BL_targ_t1(tr)) || isnan(BL_targ_t2(tr))
                error('[TF-FT][BASE] Target trial %d missing mapped baseline.', tr);
            end
            [~, lo] = min(abs(crop_times_ms - BL_targ_t1(tr)));
            [~, hi] = min(abs(crop_times_ms - BL_targ_t2(tr)));
            ii = sort([lo hi]);
            BL_idx_targ{tr} = ii(1):ii(2);
        end

       

        % ---------- Connectivity cache channels (ROI union) ----------
        freqcache_channels = unique([roi_idx.occi.all, roi_idx.pari.all, roi_idx.front.all]);
        nChanSel = numel(freqcache_channels);

        chan_hemiLR = repmat("U", nChanSel, 1);
        for k = 1:nChanSel
            ch = freqcache_channels(k);
            isL = any(ismember(ch, [roi_idx.occi.L, roi_idx.pari.L, roi_idx.front.L]));
            isR = any(ismember(ch, [roi_idx.occi.R, roi_idx.pari.R, roi_idx.front.R]));
            if isL && ~isR, chan_hemiLR(k) = "L"; end
            if isR && ~isL, chan_hemiLR(k) = "R"; end
            if isL && isR,  chan_hemiLR(k) = "U"; end
        end

        Zipsicontra = zeros(nChanSel, nTr, 'uint8');
        for tr = 1:nTr
            s = lower(strtrim(string(side_raw_targ(tr))));
            for k = 1:nChanSel
                hemi = chan_hemiLR(k);
                lab  = uint8(0);
                if hemi=="U"
                    lab = uint8(0);
                else
                    if s=="l" || s=="left" || s=="1"
                        lab = uint8(1); if hemi=="R", lab = uint8(2); end
                    elseif s=="r" || s=="right" || s=="0"
                        lab = uint8(1); if hemi=="L", lab = uint8(2); end
                    else
                        lab = uint8(0);
                    end
                end
                Zipsicontra(k,tr) = lab;
            end
        end

        % ---------- Write Zfreq cache (target) ----------
        zfreq_file_targ = fullfile(outConn, [subj_id '_targ_Zfreq.mat']);
        mf_targ = matfile(zfreq_file_targ,'Writable',true);

        Zmeta_targ = struct();
        Zmeta_targ.subject     = subj_id;
        Zmeta_targ.which       = 'targ';
        Zmeta_targ.fs          = fs;
        Zmeta_targ.times_ms    = double(times_ms);
        Zmeta_targ.freqs       = double(freqs);
        Zmeta_targ.chanlabels  = cellstr(all_labels(freqcache_channels));
        Zmeta_targ.chan_idx    = freqcache_channels;
        Zmeta_targ.roi_idx     = roi_idx;
        Zmeta_targ.chan_hemiLR = cellstr(chan_hemiLR);
        Zmeta_targ.idx_crop_ft = uint32(idx_crop_ft(:));
        Zmeta_targ.times_ms_full = double(times_ms_full(:));
        Zmeta_targ.crop_times_ms_nominal = double(crop_times_ms_nominal(:));
        mf_targ.Zmeta          = Zmeta_targ; %#ok<STRNU>

        Ztrial = struct();
        Ztrial.urevent     = uint32(urevent_targ(:));
        Ztrial.validity    = cellstr(validity_raw_targ(:));
        Ztrial.side        = cellstr(side_raw_targ(:));
        Ztrial.SOA_label   = cellstr(SOA_label_targ(:));
        Ztrial.SOA_ms      = single(SOA_raw_targ(:));
        Ztrial.RT_ms       = single(RT_ms_targ(:));
        Ztrial.acc         = single(acc_targ(:));
        Ztrial.is_uncued   = strcmpi(cellstr(validity_raw_targ(:)), 'no_cue');
        mf_targ.Ztrial     = Ztrial; %#ok<STRNU>

        % Preallocate and write ROI-union channels directly (safe dimensions)
        mf_targ.Zfreq = complex(zeros(n_f, nChanSel, nTcrop, nTr, 'single'));
        mf_targ.Zipsicontra = Zipsicontra; %#ok<STRNU>

        % Direct slice: Zfreq_crop is [freq x chan x time x trial] in EEG channel space
        % select only freqcache_channels into chanSel axis
        Ztmp = Zfreq_crop(:, freqcache_channels, :, :); % [n_f x nChanSel x nTcrop x nTr]
        mf_targ.Zfreq(1:n_f, 1:nChanSel, 1:nTcrop, 1:nTr) = Ztmp;

        % ---------- Trial-wise power accumulators ----------
        DB_band_ipsi_t   = struct();
        DB_band_contra_t = struct();
        for r = 1:numel(roi_names)
            rn = roi_names{r};
            for bb = 1:numel(band_names)
                bname = band_names{bb};
                DB_band_ipsi_t.(rn).(bname)   = zeros(nTr, nTcrop, 'single');
                DB_band_contra_t.(rn).(bname) = zeros(nTr, nTcrop, 'single');
            end
        end

        % ---------- Dense accumulators ----------
        validity_levels_t = unique(validity_raw_targ,'stable');
        side_levels_t     = ["L","R"];
        SOA_levels_t      = ["short","long","no_cue"];

        nCond = numel(validity_levels_t);
        nSide = numel(side_levels_t);
        nSOA  = numel(SOA_levels_t);
        nROI  = numel(roi_names);
        nHemi = 2;
        nCells = nCond*nSide*nSOA;

        Dense_cell_idx_t = zeros(nCond, nSide, nSOA, 'uint16');
        cnt = uint16(0);
        for ic = 1:nCond
            for isd = 1:nSide
                for isoa = 1:nSOA
                    cnt = cnt + 1;
                    Dense_cell_idx_t(ic,isd,isoa) = cnt;
                end
            end
        end

        Dense_total_sum_t = zeros(n_f, nTcrop, nROI, nHemi, nCells, 'single');
        Dense_total_n_t   = zeros(n_f, nROI, nHemi, nCells, 'uint32');
        Dense_itc_sum_t   = zeros(n_f, nTcrop, nROI, nHemi, nCells, 'single');
        Dense_itpcz_sum_t = zeros(n_f, nTcrop, nROI, nHemi, nCells, 'single');

        % ---------- Main loop over frequencies  ----------
        for fi = 1:n_f
            if fi==1 || fi==n_f || mod(fi,5)==0
                fprintf('[%s][%s] targ: fi=%d/%d, f=%.2f Hz\n', datestr(now,'HH:MM:SS'), subj_id, fi, n_f, freqs(fi));
            end

            % complex TF at this freq: [chan x time x trial]
            Zfi = squeeze(Zfreq_crop(fi, :, :, :)); % complex single [nChan x nTcrop x nTr]

            % Trial-wise dB -> band sums
            for r = 1:numel(roi_names)
                rn = roi_names{r};
                chL = roi_idx.(rn).L;
                chR = roi_idx.(rn).R;

                for tr = 1:nTr
                    if isLeft_targ(tr)
                        sel_ipsi   = chL;
                        sel_contra = chR;
                    else
                        sel_ipsi   = chR;
                        sel_contra = chL;
                    end

                    bl = BL_idx_targ{tr};

                    % IPSI
                    if ~isempty(sel_ipsi)
                        acc_db = zeros(1, nTcrop, 'single'); n_used = 0;
                        for jj = 1:numel(sel_ipsi)
                            ch = sel_ipsi(jj);
                            pw = abs(Zfi(ch,:,tr)).^2;
                            p0 = mean(pw(bl), 'omitnan');
                            if p0 <= 0 || isnan(p0), continue; end
                            db = 10*log10( pw ./ max(p0, eps) );
                            acc_db = acc_db + single(db);
                            n_used = n_used + 1;
                        end
                        if n_used > 0
                            db_ipsi = acc_db ./ n_used;
                            if any(band_sel.alpha==fi),    DB_band_ipsi_t.(rn).alpha(tr,:)    = DB_band_ipsi_t.(rn).alpha(tr,:)    + db_ipsi; end
                            if any(band_sel.beta==fi),     DB_band_ipsi_t.(rn).beta(tr,:)     = DB_band_ipsi_t.(rn).beta(tr,:)     + db_ipsi; end
                            if any(band_sel.gamma_bb==fi), DB_band_ipsi_t.(rn).gamma_bb(tr,:) = DB_band_ipsi_t.(rn).gamma_bb(tr,:) + db_ipsi; end
                        end
                    end

                    % CONTRA
                    if ~isempty(sel_contra)
                        acc_db = zeros(1, nTcrop, 'single'); n_used = 0;
                        for jj = 1:numel(sel_contra)
                            ch = sel_contra(jj);
                            pw = abs(Zfi(ch,:,tr)).^2;
                            p0 = mean(pw(bl), 'omitnan');
                            if p0 <= 0 || isnan(p0), continue; end
                            db = 10*log10( pw ./ max(p0, eps) );
                            acc_db = acc_db + single(db);
                            n_used = n_used + 1;
                        end
                        if n_used > 0
                            db_contra = acc_db ./ n_used;
                            if any(band_sel.alpha==fi),    DB_band_contra_t.(rn).alpha(tr,:)    = DB_band_contra_t.(rn).alpha(tr,:)    + db_contra; end
                            if any(band_sel.beta==fi),     DB_band_contra_t.(rn).beta(tr,:)     = DB_band_contra_t.(rn).beta(tr,:)     + db_contra; end
                            if any(band_sel.gamma_bb==fi), DB_band_contra_t.(rn).gamma_bb(tr,:) = DB_band_contra_t.(rn).gamma_bb(tr,:) + db_contra; end
                        end
                    end
                end
            end

            % Dense power + ITC/ITPCz
            PHfi = Zfi ./ max(1e-12, abs(Zfi));

            for ic = 1:nCond
                v_cur   = validity_levels_t(ic);
                tr_valid= strcmpi(validity_raw_targ, v_cur);

                for isd = 1:nSide
                    s_cur   = side_levels_t(isd);
                    tr_side = strcmpi(side_raw_targ, s_cur);

                    for isoa = 1:nSOA
                        a_cur = SOA_levels_t(isoa);

                        if  (strcmpi(v_cur,'no_cue') && ~strcmpi(a_cur,'no_cue')), continue; end
                        if (~strcmpi(v_cur,'no_cue') &&  strcmpi(a_cur,'no_cue')), continue; end

                        tr_soa = strcmpi(SOA_label_targ, a_cur);
                        tr_sel = find(tr_valid & tr_side & tr_soa);
                        if isempty(tr_sel), continue; end

                        cellIdx = Dense_cell_idx_t(ic,isd,isoa);
                        ntr = numel(tr_sel);

                        % POWER (per-frequency dB)
                        for r = 1:numel(roi_names)
                            rn  = roi_names{r};
                            chL = roi_idx.(rn).L;
                            chR = roi_idx.(rn).R;

                            % IPSI channels for this SIDE
                            if s_cur=="L", sel_ipsi = chL; else, sel_ipsi = chR; end
                            if ~isempty(sel_ipsi)
                                acc_time_sum = zeros(1, nTcrop, 'single'); nTr_used = 0;
                                for tr = tr_sel(:).'
                                    bl = BL_idx_targ{tr};
                                    acc_db = zeros(1, nTcrop, 'single'); n_used=0;
                                    for jj=1:numel(sel_ipsi)
                                        ch = sel_ipsi(jj);
                                        pw = abs(Zfi(ch,:,tr)).^2;
                                        p0 = mean(pw(bl),'omitnan');
                                        if p0<=0 || isnan(p0), continue; end
                                        acc_db = acc_db + single(10*log10(pw ./ max(p0,eps)));
                                        n_used = n_used + 1;
                                    end
                                    if n_used>0
                                        acc_time_sum = acc_time_sum + (acc_db ./ n_used);
                                        nTr_used = nTr_used + 1;
                                    end
                                end
                                if nTr_used>0
                                    Dense_total_sum_t(fi,:,r,1,cellIdx) = Dense_total_sum_t(fi,:,r,1,cellIdx) + acc_time_sum;
                                    Dense_total_n_t(  fi,  r,1,cellIdx) = Dense_total_n_t(  fi,  r,1,cellIdx) + uint32(nTr_used);
                                end
                            end

                            % CONTRA channels for this SIDE
                            if s_cur=="L", sel_contra = chR; else, sel_contra = chL; end
                            if ~isempty(sel_contra)
                                acc_time_sum = zeros(1, nTcrop, 'single'); nTr_used = 0;
                                for tr = tr_sel(:).'
                                    bl = BL_idx_targ{tr};
                                    acc_db = zeros(1, nTcrop, 'single'); n_used=0;
                                    for jj=1:numel(sel_contra)
                                        ch = sel_contra(jj);
                                        pw = abs(Zfi(ch,:,tr)).^2;
                                        p0 = mean(pw(bl),'omitnan');
                                        if p0<=0 || isnan(p0), continue; end
                                        acc_db = acc_db + single(10*log10(pw ./ max(p0,eps)));
                                        n_used = n_used + 1;
                                    end
                                    if n_used>0
                                        acc_time_sum = acc_time_sum + (acc_db ./ n_used);
                                        nTr_used = nTr_used + 1;
                                    end
                                end
                                if nTr_used>0
                                    Dense_total_sum_t(fi,:,r,2,cellIdx) = Dense_total_sum_t(fi,:,r,2,cellIdx) + acc_time_sum;
                                    Dense_total_n_t(  fi,  r,2,cellIdx) = Dense_total_n_t(  fi,  r,2,cellIdx) + uint32(nTr_used);
                                end
                            end
                        end

                        % ITC/ITPCz: per-channel across trials, then ROI mean
                        for r = 1:numel(roi_names)
                            rn  = roi_names{r};
                            chL = roi_idx.(rn).L;
                            chR = roi_idx.(rn).R;

                            for hemiTag = 1:2 % 1=ipsi, 2=contra
                                if s_cur=="L"
                                    if hemiTag==1, chan_set = chL; else, chan_set = chR; end
                                else
                                    if hemiTag==1, chan_set = chR; else, chan_set = chL; end
                                end
                                if isempty(chan_set), continue; end

                                itc_time_acc = zeros(numel(chan_set), nTcrop, 'single');
                                for jj = 1:numel(chan_set)
                                    ch = chan_set(jj);
                                    Zph = squeeze(PHfi(ch,:,tr_sel)); % [nTcrop x ntr]
                                    if size(Zph,1) ~= nTcrop
                                        Zph = reshape(Zph, [nTcrop, ntr]);
                                    end
                                    itc_time_acc(jj,:) = abs(mean(Zph, 2, 'omitnan')).';
                                end

                                itc_roi = mean(itc_time_acc, 1, 'omitnan');
                                itpcz_ch  = single(ntr) .* (itc_time_acc.^2);
                                itpcz_roi = mean(itpcz_ch, 1, 'omitnan');

                                Dense_itc_sum_t(  fi,:,r,hemiTag,cellIdx) = itc_roi;
                                Dense_itpcz_sum_t(fi,:,r,hemiTag,cellIdx) = itpcz_roi;
                            end
                        end

                    end
                end
            end
        end 

        % ---------- Finalize Dense: mean over trials ----------
        Dense_total_t = zeros(n_f, nTcrop, nROI, nHemi, nCells, 'single');
        for fi = 1:n_f
            for r = 1:nROI
                for h = 1:nHemi
                    for c = 1:nCells
                        n_here = double(Dense_total_n_t(fi,r,h,c));
                        if n_here>0
                            Dense_total_t(fi,:,r,h,c) = Dense_total_sum_t(fi,:,r,h,c) ./ n_here;
                        else
                            Dense_total_t(fi,:,r,h,c) = NaN;
                        end
                    end
                end
            end
        end

        Dense_itc_t   = Dense_itc_sum_t;
        Dense_itpcz_t = Dense_itpcz_sum_t;

        dense_file_targ = fullfile(outDense, [subj_id '_targ_dense.mat']);
        Dense = struct();
        Dense.which    = 'targ';
        Dense.freqs    = double(freqs);
        Dense.freqs_ft_raw = double(freqs_ft);
        Dense.times_ms = double(times_ms);

        Dense.levels = struct();
        Dense.levels.validity = cellstr(validity_levels_t);
        Dense.levels.side     = cellstr(side_levels_t);
        Dense.levels.SOA      = cellstr(SOA_levels_t);
        Dense.levels.ROI      = roi_names;
        Dense.levels.hemi     = {'ipsi','contra'};

        Dense.index = struct();
        Dense.index.cell = Dense_cell_idx_t;

        Dense.total = Dense_total_t;
        Dense.itc   = Dense_itc_t;
        Dense.itpcz = Dense_itpcz_t;

        save(dense_file_targ, 'Dense','-v7.3');

        % ---------- ITC rows  ----------
        itc_rows = {};
        for ic = 1:nCond
            for isd = 1:numel(side_levels_t)
                for isoa = 1:numel(SOA_levels_t)
                    cellIdx = Dense_cell_idx_t(ic,isd,isoa);
                    vlab = validity_levels_t(ic);
                    slab = side_levels_t(isd);
                    alab = SOA_levels_t(isoa);

                    if (strcmpi(vlab,'no_cue') && ~strcmpi(alab,'no_cue')), continue; end
                    if (~strcmpi(vlab,'no_cue') &&  strcmpi(alab,'no_cue')), continue; end

                    sr_norm_t = upper(strtrim(side_raw_targ));
                    tr_used = (strcmpi(validity_raw_targ,vlab) & strcmpi(sr_norm_t,slab) & strcmpi(SOA_label_targ,alab));
                    n_trials_used = sum(tr_used);
                    if n_trials_used == 0, continue; end

                    for r = 1:numel(roi_names)
                        rn = roi_names{r};
                        for hemi = ["ipsi","contra"]
                            if hemi=="ipsi", hidx = 1; else, hidx = 2; end

                            for bb = 1:numel(band_names)
                                bname = band_names{bb};
                                fidx  = band_sel.(bname);
                                if isempty(fidx), continue; end

                                itc_time   = squeeze(mean(Dense_itc_t(fidx,:,r,hidx,cellIdx), 1, 'omitnan'));
                                itpcz_time = squeeze(mean(Dense_itpcz_t(fidx,:,r,hidx,cellIdx), 1, 'omitnan'));

                                for b = 1:numel(fine_ids)
                                    idxb = fine_bins{b};
                                    itc_rows(end+1,:) = {subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemi), ...
                                        bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                        mean(itc_time(idxb),'omitnan'), mean(itpcz_time(idxb),'omitnan'), n_trials_used, numel(fidx), numel(idxb)};
                                end
                                for b = 1:numel(coll_ids)
                                    idxb = coll_bins{b};
                                    itc_rows(end+1,:) = {subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemi), ...
                                        bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                        mean(itc_time(idxb),'omitnan'), mean(itpcz_time(idxb),'omitnan'), n_trials_used, numel(fidx), numel(idxb)};
                                end
                            end
                        end
                    end
                end
            end
        end

        % ---------- Trial-wise POWER rows ----------
        trial_rows = {};
        for r = 1:numel(roi_names)
            rn = roi_names{r};
            for hemi = ["ipsi","contra"]
                for bb = 1:numel(band_names)
                    bname = band_names{bb};
                    nfb   = n_freq_in_band.(bname);

                    if hemi=="ipsi"
                        DB_here = DB_band_ipsi_t.(rn).(bname) ./ max(1,nfb);
                    else
                        DB_here = DB_band_contra_t.(rn).(bname) ./ max(1,nfb);
                    end

                    for tr = 1:nTr
                        vlab = validity_raw_targ(tr);
                        slab = side_raw_targ(tr);
                        alab = SOA_label_targ(tr);

                        if strcmp(hemi,'ipsi')
                            n_elec_nom = n_elec_hemi_nominal.(rn).ipsi.(char(slab));
                        else
                            n_elec_nom = n_elec_hemi_nominal.(rn).contra.(char(slab));
                        end

                        for b = 1:numel(fine_ids)
                            idxb = fine_bins{b};
                            trial_rows(end+1,:) = { ...
                                subj_id, EEG.setname, tr, urevent_targ(tr), "targ", char(vlab), char(slab), SOA_raw_targ(tr), char(alab), ...
                                rn, char(hemi), bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                mean(DB_here(tr, idxb), 'omitnan'), ...
                                RT_ms_targ(tr), acc_targ(tr), ...
                                n_elec_nom, nfb, numel(idxb) };
                        end

                        for b = 1:numel(coll_ids)
                            idxb = coll_bins{b};
                            trial_rows(end+1,:) = { ...
                                subj_id, EEG.setname, tr, urevent_targ(tr), "targ", char(vlab), char(slab), SOA_raw_targ(tr), char(alab), ...
                                rn, char(hemi), bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                mean(DB_here(tr, idxb), 'omitnan'), ...
                                RT_ms_targ(tr), acc_targ(tr), ...
                                n_elec_nom, nfb, numel(idxb) };
                        end
                    end
                end
            end
        end

        % =====================================================================
        % ISPC (PLV) — LONG-RANGE PHASE–PHASE CONNECTIVITY ACROSS TRIALS
        % Between-ROI electrode pairs only; within ipsi vs contra hemisphere
        % Stratified by: validity × side × SOA_label (short/long only; exclude no_cue)
        % =====================================================================
        
        % ----- Define ROI pairs -----
        roi_pairs = { ...
            'occi','pari'; ...
            'occi','front'; ...
            'pari','front' ...
        };
        
        % only compute for these validity/SOA levels
        validity_keep = ["valid","invalid"];
        SOA_keep      = ["short","long"];
        side_keep     = ["L","R"];
        
        % ISPC output rows
        ispc_rows = {};
        
        % ----- For each condition cell: compute ISPC band/timebin averages -----
        for v = 1:numel(validity_keep)
            vlab = validity_keep(v);
        
            for s = 1:numel(side_keep)
                slab = side_keep(s);
        
                for a = 1:numel(SOA_keep)
                    alab = SOA_keep(a);
        
                    % Select trials in this cell (exclude no_cue by construction)
                    tr_used = (strcmpi(validity_raw_targ, vlab) & strcmpi(side_raw_targ, slab) & strcmpi(SOA_label_targ, alab));
                    N = sum(tr_used);
        
                    if N == 0
                        continue;
                    end
       
        
                    % For each ROI pair
                    for rp = 1:size(roi_pairs,1)
                        roi1 = roi_pairs{rp,1};
                        roi2 = roi_pairs{rp,2};
        
                        % For each laterality within the selected SIDE:
                        % If slab=="L": ipsi = left hemi channels; contra = right hemi channels
                        % If slab=="R": ipsi = right hemi channels; contra = left hemi channels
                        for hemiTag = ["ipsi","contra"]
        
                            if slab=="L"
                                if hemiTag=="ipsi"
                                    chan_set_1 = roi_idx.(roi1).L;
                                    chan_set_2 = roi_idx.(roi2).L;
                                else
                                    chan_set_1 = roi_idx.(roi1).R;
                                    chan_set_2 = roi_idx.(roi2).R;
                                end
                            else
                                if hemiTag=="ipsi"
                                    chan_set_1 = roi_idx.(roi1).R;
                                    chan_set_2 = roi_idx.(roi2).R;
                                else
                                    chan_set_1 = roi_idx.(roi1).L;
                                    chan_set_2 = roi_idx.(roi2).L;
                                end
                            end
        
                            if isempty(chan_set_1) || isempty(chan_set_2)
                                continue;
                            end
        
                            n_pairs = numel(chan_set_1) * numel(chan_set_2);
                            if n_pairs == 0
                                continue;
                            end
        
                            % ----- Accumulators per band and time window -----
                            sum_ISPC  = struct();  sum_ISPCz = struct();
                            cnt_freq  = struct();
                            for bb = 1:numel(band_names)
                                bname = band_names{bb};
                                sum_ISPC.(bname)  = zeros(1, numel(fine_ids)+numel(coll_ids), 'double');
                                sum_ISPCz.(bname) = zeros(1, numel(fine_ids)+numel(coll_ids), 'double');
                                cnt_freq.(bname)  = 0;
                            end
        
                            % ----- Loop frequencies: compute pairwise ISPC(f,t) across trials, then bin -----
                            for fi = 1:n_f
        
                                % Which band?
                                in_any_band = false;
                                band_hit = "";
                                for bb = 1:numel(band_names)
                                    bname = band_names{bb};
                                    if any(band_sel.(bname) == fi)
                                        in_any_band = true;
                                        band_hit = bname;
                                        break;
                                    end
                                end
                                if ~in_any_band
                                    continue;
                                end
        
                                % Pull complex TF for selected trials at this frequency
                                chan_union = unique([chan_set_1(:); chan_set_2(:)]);
                                nSel = numel(chan_union);
                                
                                Zfi_sel = squeeze(Zfreq_crop(fi, chan_union, :, tr_used)); 
                                if ndims(Zfi_sel) ~= 3
                                    Zfi_sel = reshape(Zfi_sel, [nSel, nTcrop, N]);
                                end
                                
                                % Unit phase: PH = Z / |Z|
                                PH_f = Zfi_sel ./ max(1e-12, abs(Zfi_sel));  % [nSel x nTcrop x N]
                                
                                % Count this frequency toward the band average
                                cnt_freq.(band_hit) = cnt_freq.(band_hit) + 1;
                                
                                % Bin accumulators (sum over electrode pairs, later averaged by pair_counter)
                                nBinsTotal = numel(fine_ids) + numel(coll_ids);
                                sumPair_ISPC_bin  = zeros(1, nBinsTotal, 'double');
                                sumPair_ISPCz_bin = zeros(1, nBinsTotal, 'double');
                                
                                pair_counter = 0;
                                
                                for i1 = 1:numel(chan_set_1)
                                    chA = chan_set_1(i1);
                                    a = find(chan_union == chA, 1, 'first');
                                    if isempty(a), continue; end
                                
                                    PHa = squeeze(PH_f(a, :, :));  % [nTcrop x N]
                                    if size(PHa,1) ~= nTcrop
                                        PHa = reshape(PHa, [nTcrop, N]);
                                    end
                                
                                    for i2 = 1:numel(chan_set_2)
                                        chB = chan_set_2(i2);
                                        b = find(chan_union == chB, 1, 'first');
                                        if isempty(b), continue; end
                                
                                        PHb = squeeze(PH_f(b, :, :));  % [nTcrop x N]
                                        if size(PHb,1) ~= nTcrop
                                            PHb = reshape(PHb, [nTcrop, N]);
                                        end
                                
                                        % Phase-difference complex vector over trials, per time point
                                        D = PHa .* conj(PHb);                 % [nTcrop x N]
                                        m = mean(D, 2, 'omitnan');            % [nTcrop x 1]
                                        ispc_t  = abs(m);                     % [nTcrop x 1]
                                        ispcz_t = double(N) .* (double(ispc_t).^2);
                                
                                        % Accumulate into bins
                                        bin_idx = 0;
                                
                                        for bb2 = 1:numel(fine_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = fine_bins{bb2};
                                            if isempty(idxb), continue; end
                                            sumPair_ISPC_bin(bin_idx)  = sumPair_ISPC_bin(bin_idx)  + mean(double(ispc_t(idxb)),  'omitnan');
                                            sumPair_ISPCz_bin(bin_idx) = sumPair_ISPCz_bin(bin_idx) + mean(double(ispcz_t(idxb)), 'omitnan');
                                        end
                                
                                        for bb2 = 1:numel(coll_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = coll_bins{bb2};
                                            if isempty(idxb), continue; end
                                            sumPair_ISPC_bin(bin_idx)  = sumPair_ISPC_bin(bin_idx)  + mean(double(ispc_t(idxb)),  'omitnan');
                                            sumPair_ISPCz_bin(bin_idx) = sumPair_ISPCz_bin(bin_idx) + mean(double(ispcz_t(idxb)), 'omitnan');
                                        end
                                
                                        pair_counter = pair_counter + 1;
                                    end
                                end
        
                                if pair_counter > 0
                                    avgPair_ISPC_bin  = sumPair_ISPC_bin  ./ double(pair_counter);
                                    avgPair_ISPCz_bin = sumPair_ISPCz_bin ./ double(pair_counter);
                                else
                                    avgPair_ISPC_bin  = nan(1, nBinsTotal);
                                    avgPair_ISPCz_bin = nan(1, nBinsTotal);
                                end
        
                                sum_ISPC.(band_hit)  = sum_ISPC.(band_hit)  + avgPair_ISPC_bin;
                                sum_ISPCz.(band_hit) = sum_ISPCz.(band_hit) + avgPair_ISPCz_bin;
        
                            end 
        
                            % ----- Finalize band means and write rows -----
                            for bb = 1:numel(band_names)
                                bname = band_names{bb};
                                nfb = cnt_freq.(bname);
                                if nfb == 0
                                    continue;
                                end
        
                                ISPC_bin  = sum_ISPC.(bname)  ./ double(nfb);
                                ISPCz_bin = sum_ISPCz.(bname) ./ double(nfb);
        
                                bin_idx = 0;
        
                                for b = 1:numel(fine_ids)
                                    bin_idx = bin_idx + 1;
                                    idxb = fine_bins{b};
                                    ispc_rows(end+1,:) = { ...
                                        subj_id, "targ", roi1, roi2, char(hemiTag), char(vlab), char(slab), char(alab), ...
                                        bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                        ISPC_bin(bin_idx), ISPCz_bin(bin_idx), ...
                                        double(N), double(n_pairs), ...
                                        double(n_freq_in_band.(bname)), double(numel(idxb)) };
                                end
        
                                for b = 1:numel(coll_ids)
                                    bin_idx = bin_idx + 1;
                                    idxb = coll_bins{b};
                                    ispc_rows(end+1,:) = { ...
                                        subj_id, "targ", roi1, roi2, char(hemiTag), char(vlab), char(slab), char(alab), ...
                                        bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                        ISPC_bin(bin_idx), ISPCz_bin(bin_idx), ...
                                        double(N), double(n_pairs), ...
                                        double(n_freq_in_band.(bname)), double(numel(idxb)) };
                                end
                            end
        
                        end % hemiTag
                    end % roi pairs
                end % SOA
            end % side
        end % validity
        
        % ----- Write ISPC CSV for this subject -----
        ispc_csv = fullfile(outISPC, sprintf('%s_ispc_subject.csv', subj_id));
        if ~isempty(ispc_rows)
            writecell([ISPC_hdr; ispc_rows], ispc_csv);
            fprintf('[%s] Wrote ISPC CSV (target only): %s (rows=%d)\n', datestr(now,'HH:MM:SS'), ispc_csv, size(ispc_rows,1));
        else
            writecell(ISPC_hdr, ispc_csv);
            fprintf('[%s][WARN] No ISPC rows for %s. Created header-only CSV.\n', datestr(now,'HH:MM:SS'), subj_id);
        end

        % ---------- CSV export  ----------
        trial_csv = fullfile(outTrial, sprintf('%s_trialwise.csv', subj_id));
        if ~isempty(trial_rows)
            writecell([TW_hdr; trial_rows], trial_csv);
            fprintf('[%s] Wrote trial-wise CSV (target only): %s (rows=%d)\n', datestr(now,'HH:MM:SS'), trial_csv, size(trial_rows,1));
        else
            writecell(TW_hdr, trial_csv);
            fprintf('[%s][WARN] No trial-wise rows for %s. Created header-only CSV.\n', datestr(now,'HH:MM:SS'), subj_id);
        end

        itc_csv = fullfile(outITC, sprintf('%s_itc_subject.csv', subj_id));
        if ~isempty(itc_rows)
            writecell([ITC_hdr; itc_rows], itc_csv);
            fprintf('[%s] Wrote ITC CSV (target only): %s (rows=%d)\n', datestr(now,'HH:MM:SS'), itc_csv, size(itc_rows,1));
        else
            writecell(ITC_hdr, itc_csv);
            fprintf('[%s][WARN] No ITC rows for %s. Created header-only CSV.\n', datestr(now,'HH:MM:SS'), subj_id);
        end

        % =====================================================================
        % MATCHED-TRIAL-N EXPORTS
        % - for phase-based measures (ITPC & ISPC) - sensitive to trial
        %   number
        % - Matched within each cell: validity × side × SOA_label
        % - N* = 20 trials per cell; nIter = 100 resamples; mean across iterations
        % - Writes separate outputs with suffix: _matchedN
        % =====================================================================
        
        DO_MATCHEDN = true;   
        
        if DO_MATCHEDN
        
            fprintf('[%s][%s][MATCHEDN] Starting matched-N ITPC/ISPC block...\n', datestr(now,'HH:MM:SS'), subj_id);
        
            % -----------------------------
            % Fixed design parameters
            % -----------------------------
            N_STAR = 20;      % matched trials per cell
            N_ITER = 100;     % resampling iterations
        
            % -----------------------------
            % Reproducible RNG per subject
            % (simple deterministic seed derived from subj_id)
            % -----------------------------
            sid_chars = char(subj_id);
            sid_sum = 0;
            for kk = 1:numel(sid_chars)
                sid_sum = sid_sum + double(sid_chars(kk));
            end
            rng(100000 + sid_sum, 'twister');
        
            % -----------------------------
            % Output headers (matched-N)
            % Mirror MAIN structure but:
            % - ITPC only (no ITPCz)
            % - ISPC only (no ISPCz)
            % - always include n_trials_used (=20 when available), n_iter (=100)
            % -----------------------------
            ITC_hdr_matched = {'subj_id','which','validity','side','SOA_label','ROI','hemi','band','time_window', ...
                               'tbin_start_ms','tbin_end_ms','ITPC','n_trials_used','n_iter','n_freq_in_band','n_time_in_bin'};
        
            ISPC_hdr_matched = {'subj_id','which','roi1','roi2','laterality','validity','side','SOA_label', ...
                                'band','time_window','tbin_start_ms','tbin_end_ms', ...
                                'ISPC','n_trials_used','n_iter','n_pairs','n_freq_in_band','n_time_in_bin'};
        
            itc_rows_matched  = {};
            ispc_rows_matched = {};

            DO_DEBUG_MATCHEDN = false;
            DO_DEBUG_COMPARE_MAIN_MATCHED = false;
        
            % ================================================================
            % (A) MATCHED-N ITPC (within ROI, ipsi/contra) — same cell structure
            % Cells: validity_levels_t × side_levels_t × SOA_levels_t
            % IMPORTANT:
            % - Does NOT pool across SOA (could alter phase coherence)
            % - Does NOT pool across side 
            % - For invalid combinations (e.g., valid + no_cue SOA), skip
            %   (like pipeline above does)
            % - If N_cell < 20 => write NA rows (NaN) but keep the rows
            % ================================================================
        
            for ic = 1:numel(validity_levels_t)
                vlab = validity_levels_t(ic);
        
                for isd = 1:numel(side_levels_t)
                    slab = side_levels_t(isd);
        
                    for isoa = 1:numel(SOA_levels_t)
                        alab = SOA_levels_t(isoa);
        
                        % Mirror the SAME validity/SOA compatibility rules as main block
                        if  (strcmpi(vlab,'no_cue') && ~strcmpi(alab,'no_cue')), continue; end
                        if (~strcmpi(vlab,'no_cue') &&  strcmpi(alab,'no_cue')), continue; end
        
                        tr_cell = find( strcmpi(validity_raw_targ, vlab) & strcmpi(side_raw_targ, slab) & strcmpi(SOA_label_targ, alab) );
                        N_cell = numel(tr_cell);
        
                        % compute matched ITPC for ALL ROI × hemi × band × time bins
                        for r = 1:numel(roi_names)
                            rn = roi_names{r};
        
                            for hemiTag = ["ipsi","contra"]
        
                                % For target-locked laterality, keep exact same rule as main ITC:
                                % laterality is defined by the TARGET SIDE (slab), not by trial-wise isLeft_targ.
                                if slab=="L"
                                    if hemiTag=="ipsi"
                                        chan_set = roi_idx.(rn).L;
                                    else
                                        chan_set = roi_idx.(rn).R;
                                    end
                                else
                                    if hemiTag=="ipsi"
                                        chan_set = roi_idx.(rn).R;
                                    else
                                        chan_set = roi_idx.(rn).L;
                                    end
                                end
        
                                % If ROI has no channels for this hemi: write NA rows for ALL bands (keep full row grid)
                                if isempty(chan_set)

                                    for bb = 1:numel(band_names)
                                        bname = band_names{bb};
                                        fidx = band_sel.(bname);
                                        nfb  = double(numel(fidx));

                                        % ----- fine bins -----
                                        for b = 1:numel(fine_ids)
                                            idxb = fine_bins{b};
                                            itc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemiTag), ...
                                                bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                                NaN, double(N_STAR), double(N_ITER), nfb, double(numel(idxb)) };
                                        end

                                        % ----- collapsed bins -----
                                        for b = 1:numel(coll_ids)
                                            idxb = coll_bins{b};
                                            itc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemiTag), ...
                                                bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                                NaN, double(N_STAR), double(N_ITER), nfb, double(numel(idxb)) };
                                        end
                                    end

                                    continue;
                                end
        
                                % Loop bands (alpha/beta/gamma_bb) 
                                for bb = 1:numel(band_names)
                                    bname = band_names{bb};
                                    fidx = band_sel.(bname);
                                    nfb = numel(fidx);
        
                                    % Prepare per-bin accumulator across iterations:
                                    nBinsTotal = numel(fine_ids) + numel(coll_ids);
                                    itpc_bin_acc = zeros(1, nBinsTotal, 'double');
        
                                    if N_cell < N_STAR
                                        % Not enough trials => NA rows but KEEP structure
                                        bin_idx = 0;
        
                                        for b = 1:numel(fine_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = fine_bins{b};
                                            itc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemiTag), ...
                                                bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                                NaN, double(N_STAR), double(N_ITER), double(nfb), double(numel(idxb)) };
                                        end
                                        for b = 1:numel(coll_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = coll_bins{b};
                                            itc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemiTag), ...
                                                bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                                NaN, double(N_STAR), double(N_ITER), double(nfb), double(numel(idxb)) };
                                        end
        
                                    else
                                        % Enough trials: do resampling iterations
                                        for it = 1:N_ITER
        
                                            % Sample N_STAR trials WITHOUT replacement within this cell
                                            rp = randperm(N_cell, N_STAR);
                                            tr_sel = tr_cell(rp);

                                            % ============================================================
                                            % DEBUG: verify resampling is actually changing
                                            % and verify per-iteration ITPC values differ
                                            % Turn off by setting DO_DEBUG_MATCHEDN = false
                                            % ============================================================
                                            
                                            % Choose ONE fixed cell/ROI/hemi/band to debug
                                            DBG_validity = "valid";
                                            DBG_side     = "L";
                                            DBG_SOA      = "short";
                                            DBG_ROI      = "occi";
                                            DBG_hemi     = "ipsi";
                                            DBG_band     = "alpha";
                                            
                                            % Print only first 5 iterations
                                            if DO_DEBUG_MATCHEDN
                                                is_dbg_cell = (strcmpi(char(vlab), char(DBG_validity)) && slab==DBG_side && alab==DBG_SOA && strcmpi(rn, char(DBG_ROI)) && hemiTag==DBG_hemi && strcmpi(bname, char(DBG_band)));
                                            
                                                if is_dbg_cell
                                            
                                                    if it == 1
                                                        fprintf('[MATCHEDN-DBG][%s] DEBUG CELL: %s/%s/%s ROI=%s hemi=%s band=%s\n', subj_id, char(vlab), char(slab), char(alab), rn, char(hemiTag), bname);
                                                        fprintf('[MATCHEDN-DBG][%s] N_cell=%d, N_STAR=%d\n', subj_id, N_cell, N_STAR);
                                                        fprintf('[MATCHEDN-DBG][%s] tr_cell first 10: %s\n', subj_id, mat2str(tr_cell(1:min(10,end))'));
                                                    end
                                            
                                                    if it <= 5
                                                        fprintf('[MATCHEDN-DBG][%s] iter %d tr_sel first 10: %s\n', subj_id, it, mat2str(tr_sel(1:min(10,end))'));
                                                    end
                                            
                                                end
                                            end

                                            % Compute ITPC per frequency and then average across:
                                            % - channels (within ROI/hemi)
                                            % - freqs (within band)
                                            % - time points (within time bin)
                                            %
                                            % for each bin, accumulate a band-average ITPC value.
        
                                            % Bin-level temp accumulator for this iteration
                                            itpc_bin_thisIter = zeros(1, nBinsTotal, 'double');
        
                                            % ----- loop frequencies in this band -----
                                            for fi_i = 1:numel(fidx)
                                                fi = fidx(fi_i);
        
                                                % Pull complex TF for needed channels, time, trials:
                                                % Z: [nChanHere x nTcrop x N_STAR]
                                                Z = squeeze(Zfreq_crop(fi, chan_set, :, tr_sel));
        
                                                % Ensure 3D
                                                if ndims(Z) ~= 3
                                                    Z = reshape(Z, [numel(chan_set), nTcrop, N_STAR]);
                                                end
        
                                                % Unit-phase
                                                PH = Z ./ max(1e-12, abs(Z));  % complex unit vectors
        
                                                % ----- bins (fine + coll) -----
                                                bin_idx = 0;
        
                                                for b = 1:numel(fine_ids)
                                                    bin_idx = bin_idx + 1;
                                                    idxb = fine_bins{b};
        
                                                    if isempty(idxb)
                                                        continue;
                                                    end
        
                                                    % ITPC per channel/time: abs(mean over trials))
                                                    % PH is [nChan x nTime x N_STAR]
                                                    % compute ITPC_t as [nChan x nTimeBin]
                                                    PHbin = PH(:, idxb, :);
                                                    m = mean(PHbin, 3, 'omitnan');      % [nChan x nTimeBin]
                                                    itpc_t = abs(m);                    % [nChan x nTimeBin]
        
                                                    % average over timebin then channels
                                                    itpc_val = mean(mean(double(itpc_t), 2, 'omitnan'), 1, 'omitnan');
        
                                                    itpc_bin_thisIter(bin_idx) = itpc_bin_thisIter(bin_idx) + itpc_val;
                                                end
        
                                                for b = 1:numel(coll_ids)
                                                    bin_idx = bin_idx + 1;
                                                    idxb = coll_bins{b};
        
                                                    if isempty(idxb)
                                                        continue;
                                                    end
        
                                                    PHbin = PH(:, idxb, :);
                                                    m = mean(PHbin, 3, 'omitnan');
                                                    itpc_t = abs(m);
        
                                                    itpc_val = mean(mean(double(itpc_t), 2, 'omitnan'), 1, 'omitnan');
        
                                                    itpc_bin_thisIter(bin_idx) = itpc_bin_thisIter(bin_idx) + itpc_val;
                                                end
                                            end 
        
                                            % Average across freqs in band for this iteration
                                            if nfb > 0
                                                itpc_bin_thisIter = itpc_bin_thisIter ./ double(nfb);
                                            else
                                                itpc_bin_thisIter(:) = NaN;
                                            end

                                            if DO_DEBUG_MATCHEDN
                                                if is_dbg_cell
                                                    % report the full 0–500 ms bin (coll_id 999) if present
                                                    % bins are: fine then coll, so 999 is in the coll part.
                                                    ix999_local = find(coll_ids == 999, 1, 'first');
                                                    if ~isempty(ix999_local)
                                                        idx_999_in_total = numel(fine_ids) + ix999_local;
                                                        if it <= 5
                                                            fprintf('[MATCHEDN-DBG][%s] iter %d ITPC(0-500) = %.6f\n', subj_id, it, itpc_bin_thisIter(idx_999_in_total));
                                                        end
                                                    end
                                                end
                                            end

                                            if DO_DEBUG_MATCHEDN && is_dbg_cell && it==1
                                                % MAIN ITPC for the same cell (using ALL trials in cell)
                                                tr_all = tr_cell;         % all trials in cell
                                                N_all  = numel(tr_all);
                                            
                                                % compute ITPC(0-500) using the SAME estimator as matched (but N_all trials)
                                                ix999_local = find(coll_ids == 999, 1, 'first');
                                                idxb_999 = coll_bins{ix999_local};
                                            
                                                itpc_999_sum = 0;
                                                nfb_local = numel(fidx);
                                            
                                                for fi_i = 1:numel(fidx)
                                                    fi = fidx(fi_i);
                                            
                                                    Zfi = squeeze(Zfreq_crop(fi, :, :, :));  % [nChan x nTcrop x nTr]
                                                    PHfi = Zfi ./ max(1e-12, abs(Zfi));
                                            
                                                    % chan_set already defined for this ROI/hemi
                                                    PHc = squeeze(PHfi(chan_set, :, tr_all)); % [nChanSet x nTcrop x N_all]
                                                    m = mean(PHc, 3, 'omitnan');              % average across trials -> [nChanSet x nTcrop]
                                                    itpc_ch_t = abs(m);                        % [nChanSet x nTcrop]
                                                    itpc_roi_t = mean(itpc_ch_t, 1, 'omitnan'); % [1 x nTcrop]
                                                    itpc_999_sum = itpc_999_sum + mean(double(itpc_roi_t(idxb_999)), 'omitnan');
                                                end
                                            
                                                itpc_main_999 = itpc_999_sum / max(1, nfb_local);
                                            
                                                fprintf('[MATCHEDN-DBG][%s] MAIN (all trials) ITPC(0-500) = %.6f (N_all=%d)\n', subj_id, itpc_main_999, N_all);
                                            end
        
                                            % Accumulate across iterations
                                            itpc_bin_acc = itpc_bin_acc + itpc_bin_thisIter;
        
                                        end % iteration loop
        
                                        % Mean across iterations
                                        itpc_bin_mean = itpc_bin_acc ./ double(N_ITER);

                                        % ============================================================
                                        % DEBUG CHECK: MAIN(all trials) vs MATCHED(mean over iterations)
                                        % ============================================================
                                                                                
                                        if DO_DEBUG_COMPARE_MAIN_MATCHED
                                        
                                            % Choose one debug cell
                                            DBG_validity = "valid";
                                            DBG_side     = "L";
                                            DBG_SOA      = "short";
                                            DBG_ROI      = "occi";
                                            DBG_hemi     = "ipsi";
                                            DBG_band     = "alpha";
                                        
                                            is_dbg_cell = (strcmpi(char(vlab), char(DBG_validity)) && slab==DBG_side && alab==DBG_SOA && ...
                                                           strcmpi(rn, char(DBG_ROI)) && hemiTag==DBG_hemi && strcmpi(bname, char(DBG_band)));
                                        
                                            if is_dbg_cell
                                        
                                                % Which bin is 0–500 ms?
                                                ix999_local = find(coll_ids == 999, 1, 'first');
                                                if isempty(ix999_local)
                                                    error('[MATCHEDN-DBG] coll_ids missing 999.');
                                                end
                                                idx_999_in_total = numel(fine_ids) + ix999_local;  % bins are fine then coll
                                        
                                                % ---- Matched mean for 0–500 ----
                                                matched_999 = itpc_bin_mean(idx_999_in_total);
                                        
                                                % ---- MAIN (all trials) computed in the SAME way ----
                                                % compute the same band-averaged ITPC, but using ALL trials in tr_cell.
                                                idxb999 = coll_bins{ix999_local};   % time indices on crop grid for 0–500
                                                main_bin_acc = 0;                   % accumulate over freqs in band
                                        
                                                for fi_i = 1:numel(fidx)
                                                    fi = fidx(fi_i);
                                        
                                                    % Z: [nChanHere x nTimeBin x N_cell]
                                                    Zall = squeeze(Zfreq_crop(fi, chan_set, idxb999, tr_cell));
                                        
                                                    % Ensure 3D: [nChan x nTimeBin x N_cell]
                                                    if ndims(Zall) ~= 3
                                                        Zall = reshape(Zall, [numel(chan_set), numel(idxb999), N_cell]);
                                                    end
                                        
                                                    PHall = Zall ./ max(1e-12, abs(Zall));          % unit phase
                                                    m_all = mean(PHall, 3, 'omitnan');              % mean over trials -> [nChan x nTimeBin]
                                                    itpc_all = abs(m_all);                          % [nChan x nTimeBin]
                                        
                                                    % average over timebin then channels
                                                    itpc_val_all = mean(mean(double(itpc_all), 2, 'omitnan'), 1, 'omitnan');
                                        
                                                    main_bin_acc = main_bin_acc + itpc_val_all;
                                                end
                                        
                                                if numel(fidx) > 0
                                                    main_999 = main_bin_acc ./ double(numel(fidx));
                                                else
                                                    main_999 = NaN;
                                                end
                                        
                                                fprintf('[MATCHEDN-CHECK][%s] CELL %s/%s/%s ROI=%s hemi=%s band=%s\n', subj_id, char(vlab), char(slab), char(alab), rn, char(hemiTag), bname);
                                                fprintf('[MATCHEDN-CHECK][%s] MAIN(all N=%d) ITPC(0-500)    = %.6f\n', subj_id, N_cell, main_999);
                                                fprintf('[MATCHEDN-CHECK][%s] MATCHED(mean N*=%d, it=%d) ITPC(0-500) = %.6f\n', subj_id, N_STAR, N_ITER, matched_999);
                                                fprintf('[MATCHEDN-CHECK][%s] delta (matched-main) = %.6f\n', subj_id, matched_999 - main_999);
                                        
                                            end
                                        end
        
                                        % Write rows (fine + coll) in SAME layout as main ITC export
                                        bin_idx = 0;
        
                                        for b = 1:numel(fine_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = fine_bins{b};
                                            itc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemiTag), ...
                                                bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                                itpc_bin_mean(bin_idx), double(N_STAR), double(N_ITER), double(nfb), double(numel(idxb)) };
                                        end
        
                                        for b = 1:numel(coll_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = coll_bins{b};
                                            itc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", char(vlab), char(slab), char(alab), rn, char(hemiTag), ...
                                                bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                                itpc_bin_mean(bin_idx), double(N_STAR), double(N_ITER), double(nfb), double(numel(idxb)) };
                                        end
        
                                    end 
                                end 
                            end 
                        end 
                    end
                end 
            end 
        
            % ================================================================
            % (B) MATCHED-N ISPC (PLV) — same ROI pairs / laterality logic
            % Stratified by: valid/invalid × side × SOA(short/long)
            % IMPORTANT:
            % - Do NOT pool across SOA
            % - Do NOT pool across side
            % - Sampling performed within each (validity × side × SOA) cell
            % - If N_cell < 20 => NA rows but KEEP structure
            % ================================================================
        
            roi_pairs = { ...
                'occi','pari'; ...
                'occi','front'; ...
                'pari','front' ...
            };
        
            validity_keep = ["valid","invalid"];
            SOA_keep      = ["short","long"];
            side_keep     = ["L","R"];
        
            for v = 1:numel(validity_keep)
                vlab = validity_keep(v);
        
                for s = 1:numel(side_keep)
                    slab = side_keep(s);
        
                    for a = 1:numel(SOA_keep)
                        alab = SOA_keep(a);
        
                        tr_cell = find( strcmpi(validity_raw_targ, vlab) & strcmpi(side_raw_targ, slab) & strcmpi(SOA_label_targ, alab) );
                        N_cell = numel(tr_cell);
        
                        for rp = 1:size(roi_pairs,1)
                            roi1 = roi_pairs{rp,1};
                            roi2 = roi_pairs{rp,2};
        
                            for hemiTag = ["ipsi","contra"]
        
                                % same laterality mapping as main ISPC block (based on slab)
                                if slab=="L"
                                    if hemiTag=="ipsi"
                                        chan_set_1 = roi_idx.(roi1).L;
                                        chan_set_2 = roi_idx.(roi2).L;
                                    else
                                        chan_set_1 = roi_idx.(roi1).R;
                                        chan_set_2 = roi_idx.(roi2).R;
                                    end
                                else
                                    if hemiTag=="ipsi"
                                        chan_set_1 = roi_idx.(roi1).R;
                                        chan_set_2 = roi_idx.(roi2).R;
                                    else
                                        chan_set_1 = roi_idx.(roi1).L;
                                        chan_set_2 = roi_idx.(roi2).L;
                                    end
                                end
        
                                n_pairs = double(numel(chan_set_1) * numel(chan_set_2));
        
                                for bb = 1:numel(band_names)
                                    bname = band_names{bb};
                                    fidx = band_sel.(bname);
                                    nfb  = double(numel(fidx));
        
                                    nBinsTotal = numel(fine_ids) + numel(coll_ids);
                                    ispc_bin_acc = zeros(1, nBinsTotal, 'double');
        
                                    % If no channels or not enough trials => NA rows but keep structure
                                    if isempty(chan_set_1) || isempty(chan_set_2) || n_pairs==0 || N_cell < N_STAR
        
                                        bin_idx = 0;
                                        for b = 1:numel(fine_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = fine_bins{b};
                                            ispc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", roi1, roi2, char(hemiTag), char(vlab), char(slab), char(alab), ...
                                                bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                                NaN, double(N_STAR), double(N_ITER), double(n_pairs), double(nfb), double(numel(idxb)) };
                                        end
                                        for b = 1:numel(coll_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = coll_bins{b};
                                            ispc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", roi1, roi2, char(hemiTag), char(vlab), char(slab), char(alab), ...
                                                bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                                NaN, double(N_STAR), double(N_ITER), double(n_pairs), double(nfb), double(numel(idxb)) };
                                        end
        
                                    else
                                        % Enough trials: resample and average
                                        for it = 1:N_ITER
        
                                            rp2 = randperm(N_cell, N_STAR);
                                            tr_sel = tr_cell(rp2);
        
                                            % Bin accumulator for this iteration (band-averaged later)
                                            ispc_bin_thisIter = zeros(1, nBinsTotal, 'double');
        
                                            for fi_i = 1:numel(fidx)
                                                fi = fidx(fi_i);
        
                                                % Work in local union of channels to keep indexing safe/fast
                                                chan_union = unique([chan_set_1(:); chan_set_2(:)]);
                                                nSel = numel(chan_union);
        
                                                Z = squeeze(Zfreq_crop(fi, chan_union, :, tr_sel));  % [nSel x nTcrop x N_STAR]
                                                if ndims(Z) ~= 3
                                                    Z = reshape(Z, [nSel, nTcrop, N_STAR]);
                                                end
        
                                                PH = Z ./ max(1e-12, abs(Z));  % unit phase
        
                                                % For each pair, compute PLV across trials at each timepoint,
                                                % then average within bins; finally average over electrode pairs.
                                                sumPair_ISPC_bin = zeros(1, nBinsTotal, 'double');
                                                pair_counter = 0;
        
                                                for i1 = 1:numel(chan_set_1)
                                                    chA = chan_set_1(i1);
                                                    aidx = find(chan_union == chA, 1, 'first');
                                                    if isempty(aidx), continue; end
                                                    PHa = squeeze(PH(aidx, :, :));  % [nTcrop x N_STAR]
                                                    if size(PHa,1) ~= nTcrop
                                                        PHa = reshape(PHa, [nTcrop, N_STAR]);
                                                    end
        
                                                    for i2 = 1:numel(chan_set_2)
                                                        chB = chan_set_2(i2);
                                                        bidx = find(chan_union == chB, 1, 'first');
                                                        if isempty(bidx), continue; end
                                                        PHb = squeeze(PH(bidx, :, :)); % [nTcrop x N_STAR]
                                                        if size(PHb,1) ~= nTcrop
                                                            PHb = reshape(PHb, [nTcrop, N_STAR]);
                                                        end
        
                                                        D = PHa .* conj(PHb);           % [nTcrop x N_STAR]
                                                        m = mean(D, 2, 'omitnan');      % [nTcrop x 1]
                                                        plv_t = abs(m);                 % [nTcrop x 1]
        
                                                        bin_idx = 0;
                                                        for b = 1:numel(fine_ids)
                                                            bin_idx = bin_idx + 1;
                                                            idxb = fine_bins{b};
                                                            if isempty(idxb), continue; end
                                                            sumPair_ISPC_bin(bin_idx) = sumPair_ISPC_bin(bin_idx) + mean(double(plv_t(idxb)), 'omitnan');
                                                        end
                                                        for b = 1:numel(coll_ids)
                                                            bin_idx = bin_idx + 1;
                                                            idxb = coll_bins{b};
                                                            if isempty(idxb), continue; end
                                                            sumPair_ISPC_bin(bin_idx) = sumPair_ISPC_bin(bin_idx) + mean(double(plv_t(idxb)), 'omitnan');
                                                        end
        
                                                        pair_counter = pair_counter + 1;
                                                    end
                                                end
        
                                                if pair_counter > 0
                                                    avgPair_ISPC_bin = sumPair_ISPC_bin ./ double(pair_counter);
                                                else
                                                    avgPair_ISPC_bin = nan(1, nBinsTotal);
                                                end
        
                                                ispc_bin_thisIter = ispc_bin_thisIter + avgPair_ISPC_bin;
        
                                            end % frequency loop
        
                                            % Average across freqs in band
                                            if nfb > 0
                                                ispc_bin_thisIter = ispc_bin_thisIter ./ double(nfb);
                                            else
                                                ispc_bin_thisIter(:) = NaN;
                                            end
        
                                            ispc_bin_acc = ispc_bin_acc + ispc_bin_thisIter;
        
                                        end % iteration loop
        
                                        ispc_bin_mean = ispc_bin_acc ./ double(N_ITER);
        
                                        % Write rows (fine + coll)
                                        bin_idx = 0;
                                        for b = 1:numel(fine_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = fine_bins{b};
                                            ispc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", roi1, roi2, char(hemiTag), char(vlab), char(slab), char(alab), ...
                                                bname, fine_ids(b), fine_win(b,1), fine_win(b,2), ...
                                                ispc_bin_mean(bin_idx), double(N_STAR), double(N_ITER), double(n_pairs), double(nfb), double(numel(idxb)) };
                                        end
                                        for b = 1:numel(coll_ids)
                                            bin_idx = bin_idx + 1;
                                            idxb = coll_bins{b};
                                            ispc_rows_matched(end+1,:) = { ...
                                                subj_id, "targ", roi1, roi2, char(hemiTag), char(vlab), char(slab), char(alab), ...
                                                bname, coll_ids(b), coll_win(b,1), coll_win(b,2), ...
                                                ispc_bin_mean(bin_idx), double(N_STAR), double(N_ITER), double(n_pairs), double(nfb), double(numel(idxb)) };
                                        end
        
                                    end 
                                end 
                            end 
                        end 
                    end 
                end 
            end 
        
            % -----------------------------
            % Write matched-N CSVs (separate filenames; do NOT touch main outputs)
            % -----------------------------
            itc_csv_matched = fullfile(outITC, sprintf('%s_itc_subject_matchedN.csv', subj_id));
            if ~isempty(itc_rows_matched)
                writecell([ITC_hdr_matched; itc_rows_matched], itc_csv_matched);
                fprintf('[%s][%s][MATCHEDN] Wrote ITPC matched-N CSV: %s (rows=%d)\n', ...
                    datestr(now,'HH:MM:SS'), subj_id, itc_csv_matched, size(itc_rows_matched,1));
            else
                writecell(ITC_hdr_matched, itc_csv_matched);
                fprintf('[%s][%s][MATCHEDN][WARN] No ITPC matched-N rows. Wrote header-only CSV.\n', ...
                    datestr(now,'HH:MM:SS'), subj_id);
            end
        
            ispc_csv_matched = fullfile(outISPC, sprintf('%s_ispc_subject_matchedN.csv', subj_id));
            if ~isempty(ispc_rows_matched)
                writecell([ISPC_hdr_matched; ispc_rows_matched], ispc_csv_matched);
                fprintf('[%s][%s][MATCHEDN] Wrote ISPC matched-N CSV: %s (rows=%d)\n', ...
                    datestr(now,'HH:MM:SS'), subj_id, ispc_csv_matched, size(ispc_rows_matched,1));
            else
                writecell(ISPC_hdr_matched, ispc_csv_matched);
                fprintf('[%s][%s][MATCHEDN][WARN] No ISPC matched-N rows. Wrote header-only CSV.\n', ...
                    datestr(now,'HH:MM:SS'), subj_id);
            end
        
            fprintf('[%s][%s][MATCHEDN] Matched-N block done.\n', datestr(now,'HH:MM:SS'), subj_id);

        
        end % DO_MATCHEDN

        % ---------- Summary row ----------
        n_trl_targ = nTr;
        n_uncued_targ = sum(strcmpi(validity_raw_targ, "no_cue"));
        n_cued_targ   = sum(~strcmpi(validity_raw_targ, "no_cue") & strlength(strtrim(validity_raw_targ))>0);

        beh_mapped = sum(~isnan(urevent_targ));
        rt_nonNaN  = sum(~isnan(RT_ms_targ));
        acc_nonNaN = sum(~isnan(acc_targ));

        summary_rows(end+1,:) = {subj_id, n_trl_targ, n_cued_targ, n_uncued_targ, beh_mapped, rt_nonNaN, acc_nonNaN};

        fprintf('[%s] Subject %s done in %.1fs\n', datestr(now,'HH:MM:SS'), subj_id, toc(t_subj));

        % ---------- Clean up per subject ----------
        clear EEG ft_data freqFT F F_single Zfreq_all Zfreq_crop Zfi PHfi ...
              mf_targ Zmeta_targ Ztrial Ztmp ...
              DB_band_ipsi_t DB_band_contra_t ...
              Dense_total_sum_t Dense_total_n_t Dense_itc_sum_t Dense_itpcz_sum_t ...
              Dense_total_t Dense_itc_t Dense_itpcz_t Dense ...
              trial_rows itc_rows ispc_rows

    catch ME
        fprintf(2,'[ERR][%s] Failed: %s\n', subj_id, ME.message);
        if ~isempty(ME.stack)
            st = ME.stack(1);
            fprintf(2,'[ERR][%s] at %s (line %d)\n', subj_id, st.name, st.line);
        end
        continue;
    end

end % subject loop

%%%% ---------- Batch summary ----------
summary_csv = fullfile(outBase, 'TF_fieldtrip_batch_summary.csv');
if ~isempty(summary_rows)
    writecell([summary_hdr; summary_rows], summary_csv);
    fprintf('[TF-FT] Wrote batch summary: %s (nSubj=%d)\n', summary_csv, size(summary_rows,1));
else
    writecell(summary_hdr, summary_csv);
    fprintf('[TF-FT][WARN] No subjects processed. Wrote header-only summary.\n');
end

diary off;
fprintf('[TF-FT] DONE.\n');