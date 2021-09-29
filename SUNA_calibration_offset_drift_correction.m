function [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_calibration_offset_drift_correction(cfg,data_proc,meta_proc,data_avg,meta_avg)
%FUNCTION SUNA_calibration_offset_drift_correction
%
%  Syntax:
%    [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_calibration_offset_drift_correction(cfg,data_proc,meta_proc,data_avg,meta_avg)
%
%  Description:
%
%  Syntax:
%
%  Description:
%
%  Examples:
%
%  References:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Initialize script
dbstop if error

% Initialize correction methods
correction_methods = {'discrete' 'negative value correction'}; %'pre-deployment measurements' 'post-recovery  measurements' added below, if available.

% Initialize offset window
offsets = struct();
offsets.reference      = {};
offsets.citation       = {};
offsets.type           = {};
offsets.cum_lampon_sec = [];
offsets.NO3_uM_true    = [];
offsets.NO3_uM_SUNA    = [];
offsets.NO3_uM_TCSS_SUNA = [];
offsets.diff_TCSS      = [];

% Set filenames
cfg.path.calfile_step1 = fullfile(cfg.datadir,[cfg.project '_calibration_absolute_offsets.mat']);       % Data processed with calibration file created after recovery (i.e. usually right before the next year's deployment)
cfg.path.calfile_step2 = fullfile(cfg.datadir,[cfg.project '_calibration_reference_measurements.mat']); % Data processed with reference data collected before/after recovery and absolute offset(s) at concentration of nitrate standard(s) 
cfg.path.calfile_step3 = fullfile(cfg.datadir,[cfg.project '_calibration_discrete_samples.mat']);       % Selected discrete samples and associated metadata

%% 1 | Calibration | Calculate TCSS data from official (.cal) calibration file taken after SUNA was recovery
% Calculate data processed with calibration file created after recovery
% (i.e. usually right before the next year's deployment) 
if exist(cfg.path.calfile_step1,'file')
  load(cfg.path.calfile_step1)
else
  [cfg,data_proc,meta_proc,calfile_offsets] = SUNA_calibration_absolute_offsets(data_avg,data_proc,meta_proc,cfg);
end

%% 2 | Calibration | Calculate TCSS data using pre-deployment & post-recovery reference data
% Data processed with reference data collected before/after recovery and
% absolute offset(s) at concentration of nitrate standard(s) 
% Problem: No CTD data with pre/post deplyment measurements. This means
% cannot perform Sakamoto et al 2009 corrections, so offset not as
% precise...

if exist(cfg.path.calfile_step2,'file')
  ref_data = load(cfg.path.calfile_step2);
  % add new fields
  reference_offsets = ref_data.meta_proc.ref_offsets;
  if isfield(ref_data.data_proc,'NO3_uM_TCSS_post')
    data_proc.NO3_uM_TCSS_post = ref_data.data_proc.NO3_uM_TCSS_post;
    meta_proc.NO3_uM_TCSS_post = ref_data.meta_proc.NO3_uM_TCSS_post;
  end
  if isfield(ref_data.data_proc,'NO3_uM_TCSS_pre')
    data_proc.NO3_uM_TCSS_pre = ref_data.data_proc.NO3_uM_TCSS_pre;
    meta_proc.NO3_uM_TCSS_pre = ref_data.meta_proc.NO3_uM_TCSS_pre;
  end
else
  %[cfg,reference_offsets,post,pre] = SUNA_calibration_reference_measurements(data_avg,data_proc,meta_proc,cfg);
  [cfg,reference_offsets,data_proc,meta_proc] = SUNA_calibration_reference_measurements(data_avg,data_proc,meta_proc,cfg);
end
% Add correction methods to cell array
if exist('reference_offsets','var')
  if isfield(reference_offsets.offsets,'postdeploy')
    correction_methods = [correction_methods 'post-recovery measurements'];
  end
  if isfield(reference_offsets.offsets,'predeploy')
    correction_methods = [correction_methods 'pre-deployment measurements'];
  end
end

%% 3 | Calibration | Discrete samples
% Select discrete samples to use as calibration points.  
if exist(cfg.path.calfile_step3,'file')
  load(cfg.path.calfile_step3,'discrete')
  cfg.cal_discrete = discrete;
else
  cfg = mooring_calibration_discrete_samples(data_avg,meta_proc,cfg);
end

%% 4 | Average segment data
% initialize new variables
%       TCSS(calfile2)  TCSS(predeploy)  TCSS(postrecover) 
vars = {'NO3_uM_TCSS2' 'NO3_uM_TCSS_pre' 'NO3_uM_TCSS_post'};
vars = vars(isfield(data_proc,vars));
% Loop through data_avg and pull out average counts and nitrate
for nseg = 1:data_avg.number_of_bursts
  rng = data_proc.burst_index == nseg & data_proc.flag <= meta_proc.flag.not_evaluated;
  for nvar = 1:numel(vars)
    vname = vars{nvar};
    if nseg == 1 %initialize variable
      data_avg.(vname) = nan(data_avg.number_of_bursts,1);
      meta_avg.(vname) = meta_proc.(vname);
    end
    data_avg.(vname)(nseg) = mean(data_proc.(vname)(rng));
  end
end

%% 5 | Smooth burst averaged variables
hr_window = 35;    % 35-hour filter window for movmedian filter
idx = isfinite(data_avg.datetime);
for nv = 1:numel(vars)
  smo_var  = [vars{nv} '_smo'];
  % Initialize
  data_avg.(smo_var)      = nan(size(data_avg.(vars{nv})));
  data_avg.(smo_var)(idx) = smoothdata(data_avg.(vars{nv})(idx),'movmedian',hours(hr_window),'SamplePoints',data_avg.datetime(idx));
  % Update metadata
  meta_avg.(smo_var)  = meta_proc.(vars{nv});
  meta_avg.(smo_var).sources = vars{nv};
  meta_avg.(smo_var).method  = 'smoothdata.m, movmedian';
  meta_avg.(smo_var).window  = [num2str(hr_window) 'hours'];
  % Also smooth original field
  if ~isfield(data_proc,smo_var)
    idx2 = isfinite(data_proc.datetime);
    % Initialize
    data_proc.(smo_var)       = nan(size(data_proc.(vars{nv})));
    data_proc.(smo_var)(idx2) = smoothdata(data_proc.(vars{nv})(idx2),'movmedian',hours(hr_window),'SamplePoints',data_proc.datetime(idx2));
    % Update metadata
    meta_proc.(smo_var)  = meta_proc.(vars{nv});
    meta_proc.(smo_var).sources = vars{nv};
    meta_proc.(smo_var).method  = 'smoothdata.m, movmedian';
    meta_proc.(smo_var).window  = [num2str(hr_window) 'hours'];
  end
end

%% 6 | Plot TCSS NO3 calculated with all available calibration files
plot_NO3_and_reference_data

%% 7 | Decide which calibration file to use for processing
fprintf('--------------------------------------------------------------------------------------------------------------\n')
fprintf('Nitrate Concentrations plotted with available calibration files\n')

fprintf('--------------------------------------------------------------------------------------------------------------\n')
done_select_cal_file = 0;
while ~done_select_cal_file
  fprintf(' Which calibration file do you want to use? \n')  
  fprintf(' You will have a chance to further correct data next\n')

  for ncal = 1:numel(refs.ref_choices)
    fprintf('  <%d> %s\n',ncal,refs.ref_choices{ncal})
  end
  fprintf('  <99> STOP\n')
  chc_cal = input('  Enter choice <1/0>: ');
  if isempty(chc_cal) || chc_cal == 0
    fprintf(' incorrect choice.. try again\n');
    done_select_cal_file = 0;
  elseif chc_cal == 99
    fprintf(' entering keyboard mode.. enter "dbcont" to continue\n');
    keyboard
    done_select_cal_file = 0;
  else
    try
      TCSS_var = refs.varname{chc_cal};
      done_select_cal_file = 1;
    catch
      done_select_cal_file = 0;
    end
  end
end %% WHILE done_select_cal_file
% Initialize "corrected" data
if done_select_cal_file
  data_avg.NO3_uM_cor = data_avg.(TCSS_var);
  meta_avg.NO3_uM_cor = meta_avg.(TCSS_var);
end

%% 8 | Decide which method to use to correct the data
DONE_CORRECTING_DATA = 0;
while ~DONE_CORRECTING_DATA
  fprintf('--------------------------------------------------------------------------------------------------------------\n')
  fprintf(' Choose how to correct data: \n')
  for nc = 1:numel(correction_methods)
    % First, need to recalculate NO3 values with cur
    fprintf('  <%d> %s\n',nc,correction_methods{nc})
  end
  fprintf('  <%d> Done correcting data\n',nc+1)
  fprintf('  <99> STOP\n')
  fprintf('--------------------------------------------------------------------------------------------------------------\n')
  chc_cor_method = input(' Enter choice: ');
  switch chc_cor_method
    case 99
      keyboard
    case nc+1
      DONE_CORRECTING_DATA = 1;
    otherwise
      switch correction_methods{chc_cor_method}
        case 'discrete'
          calculate_offset_with_discrete_samples
        case 'pre-deployment measurements'
          calculate_offsets_with_reference_measurements('predeploy','NO3_uM_TCSS_pre',correction_methods{chc_cor_method})
        case 'post-recovery measurements'
          calculate_offsets_with_reference_measurements('postdeploy','NO3_uM_TCSS_post',correction_methods{chc_cor_method})
        case 'negative value correction'
          calculate_offsets_by_negative_value_adjustment
        otherwise
          fprintf('not set up yet...\n')
          keyboard
      end
  end
  % Delete corrected data plot, if necessary
  if exist('hcor','var')
    delete(hcor)
  end
  % Correct data by applying offsets
  if ~isempty(offsets.diff_TCSS)
    correct_TCSS_data_with_offsets
    % Plot corrected data
    display_text = ['NO_{3,TCSS,calibration drift corrected} ' strjoin(offsets.type,',') ];
    hcor = plot(ax,data_avg.datenum,data_avg.NO3_uM_cor,'-','Color',[0 0.4 1],'LineWidth',4,'DisplayName',display_text);
  end

end

%% 9 | Save figure and Clean up data structure (delete unused variables)
if cfg.save_figures
  standard_printfig_highrespng(fullfile(cfg.datadir,[cfg.project '_calibration_corrected_data']));
end

fields = fieldnames(data_avg);
vars(strcmp(vars,erase(TCSS_var,'_smo'))) = []; % Do not remove the original field used for corrections, becuase that may contain useful metadata
rm_fields = fields(contains(fields,vars));
for rm = 1:numel(rm_fields)
  data_avg = rmfield(data_avg,rm_fields{rm});
  meta_avg = rmfield(meta_avg,rm_fields{rm});
  if isfield(data_proc,rm_fields{rm})
    data_proc = rmfield(data_proc,rm_fields{rm});
    meta_proc = rmfield(meta_proc,rm_fields{rm});
  end
  if isfield(data_proc,[rm_fields{rm} '_rmse'])
    data_proc = rmfield(data_proc,[rm_fields{rm} '_rmse']);
    meta_proc = rmfield(meta_proc,[rm_fields{rm} '_rmse']);
  end
end

%% FUNCTION CORRECT_TCSS_DATA_WITH_OFFSETS
  function correct_TCSS_data_with_offsets
    % OFFSETS are always calculated relative to data_avg.(TCSS_var) so that
    % different corrections can be applied if necessary. 
    % Calculate corrected NO3 concentration
    % Store original
    offsets_orig = offsets;
    
    % I) Add datenum and datetime
    offsets.datetime = data_avg.datetime(nearest(data_avg.cum_lampon_sec, offsets.cum_lampon_sec));
    offsets.datenum  = data_avg.datenum(nearest(data_avg.cum_lampon_sec, offsets.cum_lampon_sec));
    
    % II) Sort offsets by cumulative lamp on time
    [~,isort] = sort(offsets.cum_lampon_sec);
    fields = fieldnames(offsets);
    for nf = 1:numel(fields)
      offsets.(fields{nf}) = offsets.(fields{nf})(isort);
    end
    
    % III) Calculate corrected NO3 concentrations
    if numel(offsets.datenum) > 1
      segments_drift_TCSS = interp1(offsets.cum_lampon_sec,offsets.diff_TCSS,data_avg.cum_lampon_sec,'linear','extrap');% data_avg: Calculate linear drift using all offsets
      fullres_drift_TCSS = interp1(offsets.cum_lampon_sec,offsets.diff_TCSS,data_proc.cum_lampon_sec,'linear','extrap');% Full resolution data: Calculate linear drift using all offsets
      correction_meaning_text = ['Nitrate data corrected with drift calculated from offsets calculated with ' strjoin(offsets.type,', ') ' measurements'];
      correction_method_text  = [TCSS_var ' + interp1(offsets.cum_lampon_sec,offsets.diff_TCSS,data_avg.cum_lampon_sec,"linear","extrap")'];
    else % Calculate single offset correction
      % data_avg: Calculate single offset correction
      segments_drift_TCSS = offsets.diff_TCSS;% Single offset 
      % Full resolution data: Calculate single offset correction
      fullres_drift_TCSS = segments_drift_TCSS; % single offset
      correction_meaning_text = 'Nitrate data corrected with single offset';
      correction_method_text  = [TCSS_var ' + offsets.diff_TCSS'];
    end
    
    % Calculated corrected burst-averaged data
    data_avg.NO3_uM_cor  = data_avg.(TCSS_var)  + segments_drift_TCSS; % data_avg
    
    % Update metadata
    meta_avg.NO3_uM_cor           = meta_avg.(TCSS_var); % Initialize
    meta_avg.NO3_uM_cor.sources   = [TCSS_var; offsets.reference]';
    meta_avg.NO3_uM_cor.reference = offsets.citation;
    meta_avg.NO3_uM_cor.method    = correction_method_text;
    meta_avg.NO3_uM_cor.meaning   = correction_meaning_text;
    meta_avg.NO3_uM_cor.offsets   = offsets;
    meta_avg.NO3_uM_cor.discrete_info = cfg.cal_discrete;
    % Calculate corrected full resolution data
    data_proc.NO3_uM_cor = data_proc.(TCSS_var) + fullres_drift_TCSS;   % full resolution
    % Update metadata
    meta_proc.NO3_uM_cor           = meta_proc.(TCSS_var); % Initialize
    meta_proc.NO3_uM_cor.sources   = [TCSS_var; offsets.reference]';
    meta_proc.NO3_uM_cor.reference = offsets.citation;
    meta_proc.NO3_uM_cor.method    = correction_method_text;
    meta_proc.NO3_uM_cor.meaning   = correction_meaning_text;
    meta_proc.NO3_uM_cor.offsets   = offsets;
    meta_proc.NO3_uM_cor.discrete_info = cfg.cal_discrete;
  end %% FUNCTION CORRECT_TCSS_DATA_WITH_OFFSETS

%% FUNCTION calculate_offset_with_discrete_samples
  function calculate_offset_with_discrete_samples
    % OFFSETS are always calculated relative to data_avg.(TCSS_var) so that
    % different corrections can be applied if necessary.
    
    % Discrete bottle samples
    % Calculate offset between TCSS SUNA data and available discrete samples
    for ndiscrete = 1:cfg.cal_discrete.number_of_samples
      n    = num2str(ndiscrete);
      scal = ['cal' n];
      if isempty(cfg.cal_discrete.(scal).info.Citation{1})
        offsets.reference = [offsets.reference; {cfg.cal_discrete.(scal).label2}];
        offsets.citation  = [offsets.citation;  {cfg.cal_discrete.(scal).label2}];
      else
        offsets.reference = [offsets.reference; cfg.cal_discrete.(scal).info.Citation(1)];
        offsets.citation  = [offsets.citation;  cfg.cal_discrete.(scal).info.Citation(1)];
      end
      % Pull out "TRUE" nitrate concentration
      bottle_NO3          = cfg.cal_discrete.(scal).(cfg.cal_discrete.(scal).cal_var);
      
      % pull out nearest SUNA value - just use single value from segment
      idx_fin  = find(isfinite(data_avg.NO3_uM_smo) & isfinite(data_avg.(TCSS_var)));
      idx_suna = nearest(data_avg.datenum(idx_fin),cfg.cal_discrete.(scal).datenum);
      idx_suna = idx_fin(idx_suna);
      
      offsets.type             = [offsets.type;     'discrete'];
      offsets.NO3_uM_true      = [offsets.NO3_uM_true; bottle_NO3];
      offsets.diff_TCSS        = [offsets.diff_TCSS       ; bottle_NO3 - data_avg.(TCSS_var)(idx_suna)];
      offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA     ; data_avg.NO3_uM_smo(idx_suna)];
      offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.(TCSS_var)(idx_suna)];
      offsets.cum_lampon_sec   = [offsets.cum_lampon_sec  ; data_avg.cum_lampon_sec(idx_suna)];
    end
  end %% FUNCTION calculate_offset_with_discrete_samples

%% FUNCTION calculate_offsets_with_reference_measurements
  function calculate_offsets_with_reference_measurements(type,field,cor_method_text)
     % OFFSETS are always calculated relative to data_avg.(TCSS_var) so that
    % different corrections can be applied if necessary.
    % Pull out specific offsets for this type
    os = reference_offsets.offsets.(type); % os just a shorthand 
    % Initialize (simply named) variables  
    os.true     = [];
    os.measured = [];
    % zero-concentration 
    if isfield(os,'baseline_uM_true')
      os.true = [os.true os.baseline_uM_true];
    end
    if isfield(os,'baseline_uM_TCSS_measured') && isfinite(os.baseline_uM_TCSS_measured)
      os.measured = [os.measured os.baseline_uM_TCSS_measured]; % Use TCSS correct data, if avaible
    elseif isfield(os,'baseline_uM_measured')
      % Assumes TCSS correction negligible...
      % i.e. DI & Nitrate standard measurements NO3(raw) == NO3(TCSS)
      os.measured = [os.measured os.baseline_uM_measured];
    end
    
    % standard-concentration 
    if isfield(os,'standard_uM_true')
      os.true = [os.true os.standard_uM_true];
    end
    if isfield(os,'standard_uM_TCSS_measured') && isfinite(os.standard_uM_TCSS_measured)
      os.measured = [os.measured os.standard_uM_TCSS_measured]; % Use TCSS correct data, if avaible
    elseif isfield(os,'standard_uM_measured')
      % Assumes TCSS correction negligible...
      % i.e. DI & Nitrate standard measurements NO3(raw) == NO3(TCSS)
      os.measured = [os.measured os.standard_uM_measured];
    end
    % Get nearest point in SUNA data
    SUNA_idx = nearest(data_avg.cum_lampon_sec,os.cum_lampon_sec);
    
    % Pull out cal files to see if they are the same
    % % Commented out March 2021 - needs more thought/testing to work.
    %     ref_calfile = meta_proc.(field).calfile;
    %     cor_calfile = meta_proc.(TCSS_var).calfile;
    %     if ~strcmp(cor_calfile,ref_calfile)
    %       % Since calibration files are different, original offset will NOT be
    %       % accurate and thus need to account for change due to calibration.
    %       for nn = 1:numel(os.true)
    %         
    %         diff_cal = data_avg.(field)(SUNA_idx) - data_avg.(TCSS_var)(SUNA_idx); % this is the difference JUST due to different calibration files
    %         % Now account for the calibration difference in the measured value
    %         os.measured(nn) = os.measured(nn) + diff_cal;
    %       end
    %     end
    
    % Calculate "TRUE" nitrate concentration at closest time
    if numel(os.measured) > 1
      % linearly interpolate if there are multiple concentrations measured
      % May 27 2021 - added extrapolation in case data_avg.(TCSS_var)(SUNA_idx) is above os.true 
      NO3_true = interp1(os.true,os.measured,data_avg.(TCSS_var)(SUNA_idx),'linear','extrap');
    else
      % Assume offset is the same over the whole concentration range
      NO3_true = (os.true - os.measured) + data_avg.(TCSS_var)(SUNA_idx);
    end

    %% Add to offsets structure
    offsets.reference        = [offsets.reference;       cor_method_text];
    offsets.type             = [offsets.type;            cor_method_text];
    offsets.citation         = [offsets.citation;        [meta_proc.(field).refdata.csvfile ' | ' num2str(meta_proc.(field).refdata.lamphrs) ' lamp hours']];
    offsets.NO3_uM_true      = [offsets.NO3_uM_true;      NO3_true];
    offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA;      data_avg.NO3_uM_smo(SUNA_idx) ];
    offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.(TCSS_var)(SUNA_idx)];
    offsets.cum_lampon_sec   = [offsets.cum_lampon_sec;   os.cum_lampon_sec];
    offsets.diff_TCSS        = [offsets.diff_TCSS;        NO3_true - data_avg.(TCSS_var)(SUNA_idx)];

  end %% FUNCTION calculate_offsets_with_reference_measurements

%% FUNCTION calculate_offsets_by_negative_value_adjustment
  function calculate_offsets_by_negative_value_adjustment
    % OFFSETS are always calculated relative to data_avg.(TCSS_var) so that
    % different corrections can be applied if necessary.
    if any(data_avg.NO3_uM_cor < 0)
      %   During the deployment at C2 in 2010–2011, the ISUS often recorded
      %   negative values and the sensor failed in July 2011, eliminating the
      %   possibility of an in-situ calibration from the recovery CTD cast. In
      %   this instance, the calibration drift correction used the initial
      %   in-situ calibration point and the most negative daily mean value
      %   (observed on 12 June), which was set to zero. The resulting pattern was
      %   similar to the nitrate time series at C1, which had a double maximum
      %   between late March and late May, and ~9 ?M drop in nitrate on 4–7 June
      %   (not shown). After calibration, the C2 time series 2016–2017 had
      %   periods of negative values. For this time series, a secondary drift
      %   correction was applied by setting the most negative daily mean value
      %   (observed on 14 November) to zero. [Mordy et al., 2020]
      [~,imin_TCSS] = nanmin(data_avg.(TCSS_var));
      offsets.reference        = [offsets.reference;        [datestr(data_avg.datenum(imin_TCSS) ) ' Most negative corrected ' TCSS_var ' value']];
      offsets.type             = [offsets.type;             'zero'];
      offsets.citation         = [offsets.citation;         'zero'];
      offsets.cum_lampon_sec   = [offsets.cum_lampon_sec;   data_avg.cum_lampon_sec(imin_TCSS)];
      offsets.NO3_uM_true      = [offsets.NO3_uM_true;      0.0];
      offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA;      data_avg.NO3_uM_smo(imin_TCSS)];
      offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.(TCSS_var)(imin_TCSS)];
      
      % now calculate differences
      offsets.diff_TCSS = [offsets.diff_TCSS; 0.0 - data_avg.(TCSS_var)(imin_TCSS)];
    else
      fprintf('No negative values to adjust... skipping this\n')
    end
  end %% FUNCTION calculate_offsets_by_negative_value_adjustment

%% FUNCTION PLOT_NO3_AND_REFERENCE_DATA
  function plot_NO3_and_reference_data
    % add all available data to "refs" structure so can calculate offset
    % later, if chosen. 
    makefig; ax=subplot(6,1,1:4); hold(ax,'on'); grid(ax,'on'); ax.YDir = 'normal';  ax.Color = [0.75 0.75 .75]; ax.Box = 'on';
    ax.Title.String = ['SUNA ' cfg.mooring.SN.SUNA ' Calibration corrected data'];
    ax.YLabel.String = 'NO_3 [\muM]';
    axis(ax,'tight');
    ax2=subplot(6,1,5:6); hold(ax2,'on');grid(ax2,'on');ax2.YDir = 'normal'; ax2.Color = [0.75 0.75 .75]; ax2.Box = 'on';
    ax2.YLabel.String = 'Reference spectra';
    ax2.XLabel.String = 'Wavelength [nm]';
    axis(ax2,'tight');
    ax2.XLim = [190 400];
    % Plot raw and TCSS corrected data with original calibration file
    calstr1 =  ['Cal1: ' calfile_offsets.calfiles{1} ' | ' calfile_offsets.caldates{1} ' | ' num2str(calfile_offsets.lamphour(1),'%.0f') ' lamp hours'];
    refs = struct();
    refs.ref_choices = {calstr1};
    refs.varname     = {'NO3_uM_TCSS_smo'};
    refs.lamphours   = calfile_offsets.lamphour(1);
    refs.lampsec     = calfile_offsets.lampsec(1);
    plot(ax,data_avg.datenum,data_avg.NO3_uM_smo,     'k-','LineWidth',2,'DisplayName',['NO_{3(raw,cal1,smoothed) }  ' calstr1]);
    plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_smo,'g-','LineWidth',2,'DisplayName',['NO_{3(TCSS,cal1,smoothed) } ' calstr1]);
    % Plot the original calibration reference spectra
    plot(ax2,calfile_offsets.cal1.Wavelength,calfile_offsets.cal1.Reference, 'k-','LineWidth',2,'DisplayName',calstr1);
    % Plot TCSS corrected data with calibration file from next deployment
    if isfield(data_avg,'NO3_uM_TCSS2_smo')
      calstr2 =  ['Cal2: ' calfile_offsets.calfiles{2} ' | ' calfile_offsets.caldates{2} ' | ' num2str(calfile_offsets.lamphour(2),'%.0f') ' lamp hours'];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS2_smo,                      'b-','LineWidth',2,'DisplayName',['NO_{3(TCSS,cal2,smoothed) } ' calstr2]);
      plot(ax2,calfile_offsets.cal1.Wavelength,calfile_offsets.cal2.Reference, 'b-','LineWidth',2,'DisplayName',calstr2);
      %Add to refs structure
      refs.ref_choices = [refs.ref_choices; calstr2];
      refs.varname     = [refs.varname; 'NO3_uM_TCSS2_smo'];
      refs.lamphours   = [refs.lamphours; calfile_offsets.lamphour(2)];
      refs.lampsec     = [refs.lampsec;   calfile_offsets.lampsec(2)];
    end
    % Plot TCSS corrected data with pre-deployment reference spectra
    if isfield(data_avg,'NO3_uM_TCSS_pre_smo')
      calstr_pre = ['CAL_pre: ' meta_proc.NO3_uM_TCSS_pre.calfile ' | ' num2str(meta_proc.NO3_uM_TCSS_pre.refdata.lamphrs) ' lamp hours'];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_pre_smo,      'r-','LineWidth',2,'DisplayName',['NO_{3(TCSS,pre,smoothed) } ' strrep(calstr_pre,'_','\_')]);
      plot(ax2,calfile_offsets.cal1.Wavelength,meta_proc.NO3_uM_TCSS_pre.caldata.Reference, 'r-','LineWidth',2,'DisplayName', strrep(calstr_pre,'_','\_'));
      % Add to refs structure
      refs.ref_choices = [refs.ref_choices; calstr_pre];
      refs.varname     = [refs.varname; 'NO3_uM_TCSS_pre_smo'];
      refs.lamphours   = [refs.lamphours; meta_proc.NO3_uM_TCSS_pre.refdata.lamphrs];
      refs.lampsec     = [refs.lampsec;   meta_proc.NO3_uM_TCSS_pre.refdata.lampsec];
    end
    % Plot TCSS corrected data with post-recovery reference spectra
    if isfield(data_avg,'NO3_uM_TCSS_post_smo')
      calstr_post = ['CAL_post: ' meta_proc.NO3_uM_TCSS_post.calfile ' | ' num2str(meta_proc.NO3_uM_TCSS_post.refdata.lamphrs) ' lamp hours'];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_post_smo,      'm-','LineWidth',2,'DisplayName',['NO_{3(TCSS,post,smoothed) } ' strrep(calstr_post,'_','\_')]);
      plot(ax2,calfile_offsets.cal1.Wavelength,meta_proc.NO3_uM_TCSS_post.caldata.Reference, 'm-','LineWidth',2,'DisplayName', strrep(calstr_post,'_','\_'));
      % Add to refs structure
      refs.ref_choices = [refs.ref_choices; calstr_post];
      refs.varname     = [refs.varname; 'NO3_uM_TCSS_post_smo'];
      refs.lamphours   = [refs.lamphours; meta_proc.NO3_uM_TCSS_post.refdata.lamphrs];
      refs.lampsec     = [refs.lampsec;   meta_proc.NO3_uM_TCSS_post.refdata.lampsec];
    end
    datetick(ax,'x');
    hl = legend(ax,'show');  hl.FontSize = 12;
    hl.Units = 'normalized'; hl.Position(1:2) = [0.002 0.85];
    hl.Location = 'south';
    try % Issues with legend in R2019a on mac...
      axes(ax2);
      hl2 = legend(ax2,'show'); 
      hl2.Location = 'northeast';
      hl2.FontSize = 12;
    catch
      fprintf('legend not working..\n')
      keyboard
    end
    
    
    % plot discrete sample
    if isfield(cfg,'cal_discrete')
      clrs = hsv(cfg.cal_discrete.number_of_samples);
      for ndiscrete = 1:cfg.cal_discrete.number_of_samples
        n = num2str(ndiscrete); scal = ['cal' n];
        bottle_NO3 = cfg.cal_discrete.(scal).(cfg.cal_discrete.(scal).cal_var);
        plot(ax,cfg.cal_discrete.(scal).datenum,bottle_NO3,'kh','MarkerSize',20,'LineWidth',0.5,'MarkerFaceColor',clrs(ndiscrete,:),'DisplayName',cfg.cal_discrete.(scal).header)
      end
    end
    fprintf('\n');
  end %% FUNCTION PLOT_NO3_AND_REFERENCE_DATA


end %% MAIN FUNCTION