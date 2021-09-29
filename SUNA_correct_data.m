function [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_correct_data(cfg,data_proc,meta_proc,data_avg,meta_avg)
%FUNCTION SUNA_correct_data
%
%  Syntax:
%    [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_correct_data(cfg,data_proc,meta_proc,data_avg,meta_avg)
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
%% Initialize script
dbstop if error
% Initialize offset window
offsets = struct();
offsets.reference      = {};
offsets.citation       = {};
offsets.type           = {};
offsets.cum_lampon_sec = [];
offsets.NO3_uM_true    = [];
offsets.NO3_uM_SUNA    = [];
offsets.NO3_uM_TCSS_SUNA = [];
offsets.diff_raw       = [];
offsets.diff_TCSS      = [];
% initialize temporary array to store corrected data
NO3_cor                = nan(size(data_avg.datenum));

%% 1 | Calibration | previous .cal files
% Problem: No CTD data with pre/post deplyment measurements. This means
% cannot perform Sakamoto et al 2009 corrections, so offset not as
% precise...
[cfg,data_proc,meta_proc,calfile_offsets] = SUNA_calibration_absolute_offsets(data_avg,data_proc,meta_proc,cfg);

%% 2 | Calibration | pre-deployment & post-recovery reference data
% Problem: No CTD data with pre/post deplyment measurements. This means
% cannot perform Sakamoto et al 2009 corrections, so offset not as
% precise...
[cfg,reference_offsets,post,pre] = SUNA_calibration_reference_measurements(data_avg,data_proc,meta_proc,cfg);

%% 3 | Calibration | Discrete samples
[cfg,discrete_offsets] = SUNA_calibration_discrete_samples(data_avg,meta_proc,cfg);

%% 4 | Calculate TCSS Nitrate using reference measurements
if ~isempty(pre)
  % before window cleaning
  cal_pre = table2struct(calfile_offsets.cal1,'ToScalar',1); % all variables the same except will change Reference
  cal_pre.Reference = nanmean(pre.spectrum_channels(pre.standard_idx,:),1)';
  cal_pre.csvfile = pre.filename(pre.standard_idx(1),:);
  cal_pre.lamphrs = round(nanmean(pre.cum_lampon_sec(pre.standard_idx))./60./60);
  cal_pre.lampsec = nanmean(pre.cum_lampon_sec(pre.standard_idx));
  fprintf('NEED TO FIGURE OUT IF cal.NO3 IS ABSORPTION OR WHAT...\n')
  %Recalculate data with post recovery reference measurements
  fprintf('Calculating TCSS Nitrate concentration using pre-deployment reference measurements\n');
  [NTR_pre,~] = ISUS_REPROCESSOR_v2_bi(data_proc,meta_proc,cal_pre,1);
  %[NTR_pre,~] = ISUS_REPROCESSOR_v2_bi(data_proc,[],cal_pre,1);
  data_proc.NO3_uM_TCSS_pre = NTR_pre;
  fprintf(' WRITE CAL FILE\n')
  keyboard 
end

if ~isempty(post)
  cal_post = table2struct(calfile_offsets.cal1,'ToScalar',1); % all variables the same except will change Reference
  cal_post.Reference  = nanmean(post.spectrum_channels(post.baseline_idx,:),1)';
  cal_post.csvfile    = post.filename(post.baseline_idx(1),:);
  cal_post.lamphrs    = round(nanmean(post.cum_lampon_sec(post.baseline_idx))./60./60);
  cal_post.lampsec    = nanmean(post.cum_lampon_sec(post.baseline_idx));
  fprintf('NEED TO FIGURE OUT IF cal.NO3 IS ABSORPTION OR WHAT...\n')
  %Recalculate data with post recovery reference measurements
  fprintf('Calculating TCSS Nitrate concentration using post-recovery reference measurements\n');
  [NTR_post,~] = ISUS_REPROCESSOR_v2_bi(data_proc,meta_proc,cal_post,1);
  %[NTR_post,~] = ISUS_REPROCESSOR_v2_bi(data_proc,[],cal_post,1);
  data_proc.NO3_uM_TCSS_post = NTR_post;
  fprintf(' WRITE CAL FILE\n')
  keyboard
end

%% 5 | Average segment data
% initialize new variables
%       TCSS(calfile2)  TCSS(predeploy)  TCSS(postrecover) 
vars = {'NO3_uM_TCSS2' 'NO3_uM_TCSS_pre' 'NO3_uM_TCSS_post'};
vars = vars(isfield(data_proc,vars));
% Loop through data_avg and pull out average counts and nitrate
for nseg = 1:data_avg.number_of_segments
  rng = data_proc.burst_index == nseg & data_proc.flag <= meta_proc.flag.not_evaluated;
  for nvar = 1:numel(vars)
    vname = vars{nvar};
    if nseg == 1 %initialize variable
      data_avg.(vname) = nan(data_avg.number_of_segments,1);
    end
    data_avg.(vname)(nseg) = mean(data_proc.(vname)(rng));
  end
end

%% 6 | smooth new variables
hr_window = 35;    % 35-hour filter window for movmedian filter
idx = isfinite(data_avg.datetime);
for nv = 1:numel(vars)
  smo_var  = [vars{nv} '_smo'];
  % Initialize
  data_avg.(smo_var)      = nan(size(data_avg.(vars{nv})));
  data_avg.(smo_var)(idx) = smoothdata(data_avg.(vars{nv})(idx),'movmedian',hours(hr_window),'SamplePoints',data_avg.datetime(idx));
end

%% PLOT TCSS NO3 with different reference data
plot_NO3_and_reference_data

%% Decide if want to use reference measurements to correct drift
fprintf('--------------------------------------------------------------------------------------------------------------\n')
fprintf('Nitrate Concentrations plotted with available calibration files and pre/post deployment reference measurements\n')
fprintf('--------------------------------------------------------------------------------------------------------------\n')
fprintf('POST-CRUISE CAL OFFSET: Do you want to calculate offset based on one of these?\n')
chc_ref = input('  Enter choice <1/0>: ');
if isempty(chc_ref); chc_ref = 0; end % Default to NO
if chc_ref
  fprintf('\n Select which dataset to use to calculate offset\n')
  for nref = 1:numel(refs.ref_choices)
    fprintf('  <%d> %s\n',nref,refs.ref_choices{nref})
  end
  num_ref = input('  Enter choice: ');
  idx_first = find(isfinite(data_avg.NO3_uM),1);
  idx_end   = find(isfinite(data_avg.NO3_uM),1,'last');
  % add information to offsets structure
  post_offset = data_avg.(refs.varname{num_ref})(idx_end) - data_avg.NO3_uM_TCSS_smo(idx_end);
  offsets.reference        = [offsets.reference;      'post-recovery reference'];
  offsets.type             = [offsets.type;           'post-recovery reference'];
  offsets.citation         = [offsets.citation;       refs.ref_choices{nref}];
  offsets.cum_lampon_sec   = [offsets.cum_lampon_sec; refs.lampsec(nref)];
  offsets.NO3_uM_true      = [offsets.NO3_uM_true;      data_avg.(refs.varname{num_ref})(idx_end) ];
  offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA;      data_avg.NO3_uM_smo(idx_end) ];
  offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.NO3_uM_TCSS_smo(idx_end) ];
  offsets.diff_TCSS        = [offsets.diff_TCSS;        data_avg.(refs.varname{num_ref})(idx_end) - data_avg.NO3_uM_TCSS_smo(idx_end)];
  offsets.diff_raw         = [offsets.diff_raw;         data_avg.(refs.varname{num_ref})(idx_end) - data_avg.NO3_uM_smo(idx_end)];
  % calculate drift based on post-recovery reference measurements
  postcal_offset         = interp1([data_avg.cum_lampon_sec(idx_first) refs.lampsec(num_ref)],[0 post_offset],data_avg.cum_lampon_sec,'linear','extrap');
  NO3_cor                = data_avg.NO3_uM_TCSS_smo + postcal_offset;
  plot(ax,data_avg.datenum,NO3_cor, 'k--','LineWidth',3,'DisplayName','NO_{3(TCSS,corrected 1) }');

end

%% Add offset for pre deployment reference measurements
if exist('pre','var') && ~isempty(pre)
  % Assume calibration file correct... 
  % Assume TCSS correction negligible... 
  %i.e. pre-deployment DI & Nitrate standard measurements NO3(raw) == NO3(TCSS)
  idx_first = find(isfinite(data_avg.NO3_uM_TCSS_smo),1);
  NO3_TCSS_pre = interp1([0 pre.NO3_uM_standard(1)],[pre.baseline_NO3 pre.standard_NO3],data_avg.NO3_uM_TCSS_smo(idx_first));
  NO3_raw_pre  = interp1([0 pre.NO3_uM_standard(1)],[pre.baseline_NO3 pre.standard_NO3],data_avg.NO3_uM_smo(idx_first));
  diff_pre     = NO3_TCSS_pre - data_avg.NO3_uM_TCSS_smo(idx_first);
  diff_raw     = NO3_raw_pre  - data_avg.NO3_uM_smo(idx_first);
  
  offsets.reference        = [offsets.reference;      'pre-deployment reference'];
  offsets.type             = [offsets.type;           'pre-deployment reference'];
  offsets.citation         = [offsets.citation;       [cal_pre.csvfile ' | ' num2str(cal_pre.lamphrs) ' lamp hours']];
  offsets.NO3_uM_true      = [offsets.NO3_uM_true;      NO3_TCSS_pre];
  offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA;      data_avg.NO3_uM_smo(idx_first) ];
  offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.NO3_uM_TCSS_smo(idx_first) ];
  offsets.cum_lampon_sec   = [offsets.cum_lampon_sec; cal_pre.lampsec];
  offsets.diff_TCSS        = [offsets.diff_TCSS;      diff_pre];
  offsets.diff_raw         = [offsets.diff_raw;       diff_raw];
  % calculate drift based on post-recovery reference measurements
  pre_offsets            = interp1([data_avg.cum_lampon_sec(idx_first)  data_avg.cum_lampon_sec(idx_end)],[diff_pre 0],data_avg.cum_lampon_sec,'linear','extrap');
  NO3_cor_post           = NO3_cor;
  NO3_cor                = NO3_cor + pre_offsets;
  plot(ax,data_avg.datenum,NO3_cor, 'r--','LineWidth',3,'DisplayName','NO_{3(TCSS,corrected 2) }');
end


%% Discrete bottle samples
clrs = hsv(discrete_offsets.number_of_samples);
for ndiscrete = 1:discrete_offsets.number_of_samples
  n = num2str(ndiscrete);
  %offsets.datenum(ncal)  = datenum(discrete_offsets.(['Bottle_date' n]));
  %offsets.datetime(ncal) = datetime(offsets.datenum(ncal),'ConvertFrom','datenum');
  %offsets.reference{ncal} = ['Cruise:' char(discrete_offsets.(['info' n]).CRUISE{1}) ' Station:' char(discrete_offsets.(['info' n]).StationID{1}) ...
  %  ' Cast:' num2str(discrete_offsets.(['info' n]).CastID(1)) ' Niskin:' num2str(discrete_offsets.(['info' n]).Niskin(1))];
  discrete_string   =  ['Cruise:' char(discrete_offsets.(['info' n]).CRUISE{1}) ' Station:' char(discrete_offsets.(['info' n]).StationID{1}) ...
                        ' Cast:' num2str(discrete_offsets.(['info' n]).CastID(1)) ' Niskin:' num2str(discrete_offsets.(['info' n]).Niskin(1))];
  offsets.reference = [offsets.reference; discrete_string];
  offsets.citation  = [offsets.citation;  char(discrete_offsets.(['info' n]).Citation{1})];
  offsets.type      = [offsets.type;     'discrete'];
  % Use interpolated values to calculate offsets
  if any(contains(fieldnames(discrete_offsets),'interpolated'))
    try
      bottle_NO3 = discrete_offsets.(['Bottle_NO3_uM_interpolated' n]);
    catch
      bottle_NO3 = discrete_offsets.(['Bottle_NO3_NO2_uM_interpolated' n]);
    end
  elseif any(contains(fieldnames(discrete_offsets),'averaged'))
    try
      bottle_NO3 = discrete_offsets.(['Bottle_NO3_uM_averaged' n]);
    catch
      bottle_NO3 = discrete_offsets.(['Bottle_NO3_NO2_uM_averaged' n]);
    end
  else
     try
       bottle_NO3 = discrete_offsets.(['Bottle_NO3_uM' n]);
     catch
       bottle_NO3 = discrete_offsets.(['Bottle_NO3_NO2_uM' n]);
     end
  end
  offsets.NO3_uM_true = [offsets.NO3_uM_true; bottle_NO3];
  % pull out nearest SUNA value - just use single value from segment
  idx_fin  = find(isfinite(data_avg.NO3_uM_smo) & isfinite(data_avg.NO3_uM_TCSS_smo));
  idx_suna = nearest(data_avg.datenum(idx_fin),discrete_offsets.(['Bottle_datenum' n]));
  idx_suna = idx_fin(idx_suna);
  offsets.diff_raw         = [offsets.diff_raw        ; bottle_NO3 - data_avg.NO3_uM_smo(idx_suna)];
  offsets.diff_TCSS        = [offsets.diff_TCSS       ; bottle_NO3 - data_avg.NO3_uM_TCSS_smo(idx_suna)];
  offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA     ; data_avg.NO3_uM_smo(idx_suna)];
  offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.NO3_uM_TCSS_smo(idx_suna)];
  offsets.cum_lampon_sec   = [offsets.cum_lampon_sec  ; data_avg.cum_lampon_sec(idx_suna)];
  % plot discrete sample
  % calculate drift based on post-recovery reference measurements
  [~,inot_near]  = max(abs(data_avg.datenum - discrete_offsets.(['Bottle_datenum' n])));
  diff_bot       = bottle_NO3 - NO3_cor(idx_suna);
  bot_offset     = interp1([data_avg.cum_lampon_sec(inot_near)  data_avg.cum_lampon_sec(idx_suna)],[0 diff_bot],data_avg.cum_lampon_sec,'linear','extrap');
  NO3_cor_old    = NO3_cor;
  NO3_cor        = NO3_cor + bot_offset;
  % Plot data corrected with discrete bottle and discrete bottle 
  plot(ax,data_avg.datenum,NO3_cor,':','LineWidth',2,'Color',clrs(ndiscrete,:),'DisplayName',discrete_string)
  plot(ax,discrete_offsets.(['Bottle_datenum' n]),bottle_NO3,'kh','MarkerSize',10,'LineWidth',0.2,'MarkerFaceColor',clrs(ndiscrete,:),'DisplayName',discrete_string)
end

%% 4 | Negative values offset
if any(NO3_cor < 0)
  fprintf('Negative values!\n')
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
  [~,imin_TCSS] = nanmin(NO3_cor);
  offsets.reference        = [offsets.reference;        [datestr(data_avg.datenum(imin_TCSS) ) ' Most negative corrected NO3_uM_TCSS value']];
  offsets.type             = [offsets.type;             'zero'];
  offsets.citation         = [offsets.citation;         'zero'];
  offsets.cum_lampon_sec   = [offsets.cum_lampon_sec;   data_avg.cum_lampon_sec(imin_TCSS)];
  offsets.NO3_uM_true      = [offsets.NO3_uM_true;      0.0];
  offsets.NO3_uM_SUNA      = [offsets.NO3_uM_SUNA;      data_avg.NO3_uM_smo(imin_TCSS)];
  offsets.NO3_uM_TCSS_SUNA = [offsets.NO3_uM_TCSS_SUNA; data_avg.NO3_uM_TCSS_smo(imin_TCSS)];
 
  % now calculate differences
  offsets.diff_raw  = [offsets.diff_raw;  0.0  - data_avg.NO3_uM_smo(imin_TCSS)];
  offsets.diff_TCSS = [offsets.diff_TCSS; 0.0 - data_avg.NO3_uM_TCSS_smo(imin_TCSS)];
end

%% 5 | Calculate corrected NO3 concentration
% first sort offsets by cumulative lamp on time
offsets2 = offsets;
[~,isort] = sort(offsets.cum_lampon_sec);
fields = fieldnames(offsets);
for nf = 1:numel(fields)
  offsets.(fields{nf}) = offsets.(fields{nf})(isort);
end
% add datenum and datetime
idx_seg = nearest(data_avg.cum_lampon_sec, offsets.cum_lampon_sec);
offsets.datetime = data_avg.datetime(nearest(data_avg.cum_lampon_sec, offsets.cum_lampon_sec));
offsets.datenum  = data_avg.datenum(nearest(data_avg.cum_lampon_sec, offsets.cum_lampon_sec));

%% 6 | Calculate corrected NO3 concentrations
if numel(offsets.datenum) > 1
  % data_avg: Calculate linear drift using all offsets
  segments_drift_raw =  interp1(offsets.cum_lampon_sec,offsets.diff_raw, data_avg.cum_lampon_sec,'linear','extrap');
  segments_drift_TCSS = interp1(offsets.cum_lampon_sec,offsets.diff_TCSS,data_avg.cum_lampon_sec,'linear','extrap');
  % Full resolution data: Calculate linear drift using all offsets
  fullres_drift_raw =  interp1(offsets.cum_lampon_sec,offsets.diff_raw, data_proc.cum_lampon_sec,'linear','extrap');
  fullres_drift_TCSS = interp1(offsets.cum_lampon_sec,offsets.diff_TCSS,data_proc.cum_lampon_sec,'linear','extrap');
  correction_meaning_text = ['Nitrate data corrected with drift calculated from offsets calculated with ' strjoin(offsets.type,', ') ' measurements'];
  correction_method_text  = 'NO3_uM_TCSS_smo + interp1(offsets.cum_lampon_sec,offsets.diff_TCSS,data_avg.cum_lampon_sec,"linear","extrap")';
else % Calculate single offset correction
  % data_avg: Calculate single offset correction
  segments_drift_raw =  offsets.diff_raw; % Single offset
  segments_drift_TCSS = offsets.diff_TCSS;% Single offset
  % Full resolution data: Calculate single offset correction
  fullres_drift_raw  = segments_drift_raw;  % single offset
  fullres_drift_TCSS = segments_drift_TCSS; % single offset
  correction_meaning_text = 'Nitrate data corrected with single offset';
  correction_method_text  = 'NO3_uM_TCSS_smo + offsets.diff_TCSS';
end

% Calculated corrected data
data_avg.NO3_uM_cor  = data_avg.NO3_uM_TCSS_smo  + segments_drift_TCSS; % data_avg
data_proc.NO3_uM_cor = data_proc.NO3_uM_TCSS_smo + fullres_drift_TCSS;   % full resolution
% Update metadata  
meta_proc.NO3_uM_cor           = meta_proc.NO3_uM_TCSS; % Initialize 
meta_proc.NO3_uM_cor.sources   = ['NO3_uM_TCSS_smo'; offsets.reference]';
meta_proc.NO3_uM_cor.reference = offsets.citation;
meta_proc.NO3_uM_cor.method    = correction_method_text;
meta_proc.NO3_uM_cor.meaning   = correction_meaning_text;
meta_proc.NO3_uM_cor.offsets   = offsets;

%% Plot single axes
% first save larger figure
if cfg.save_figures
  savename = fullfile(cfg.datadir,[cfg.project '_correct_data_options1']);
  standard_printfig_highrespng(savename);
end
plot_NO3
plot(ax,data_avg.datenum,data_avg.NO3_uM_cor,'-','Color',[0 0.4 1],'LineWidth',4,'DisplayName','NO_{3,TCSS,calibration drift corrected}')
% Plot bottle data
for ndiscrete = 1:discrete_offsets.number_of_samples
  n = num2str(ndiscrete);
  discrete_string   =  ['Cruise:' char(discrete_offsets.(['info' n]).CRUISE{1}) ' Station:' char(discrete_offsets.(['info' n]).StationID{1}) ...
                        ' Cast:' num2str(discrete_offsets.(['info' n]).CastID(1)) ' Niskin:' num2str(discrete_offsets.(['info' n]).Niskin(1))];
  % Use interpolated values to calculate offsets
  if any(contains(fieldnames(discrete_offsets),'interpolated'))
    try bottle_NO3 = discrete_offsets.(['Bottle_NO3_uM_interpolated' n]);
    catch; bottle_NO3 = discrete_offsets.(['Bottle_NO3_NO2_uM_interpolated' n]);
    end
  else
     try bottle_NO3 = discrete_offsets.(['Bottle_NO3_uM' n]);
     catch; bottle_NO3 = discrete_offsets.(['Bottle_NO3_NO2_uM' n]);
     end
  end
  plot(ax,discrete_offsets.(['Bottle_datenum' n]),bottle_NO3,'kh','MarkerSize',12,'LineWidth',0.2,'MarkerFaceColor',clrs(ndiscrete,:),'DisplayName',discrete_string)
end
if cfg.save_figures
  savename = fullfile(cfg.datadir,[cfg.project '_offset_drift_corrected_data']);
  standard_printfig_highrespng(savename);
end
%% Clean up data structure - delete temporary variables
fields = fieldnames(data_avg);
rm_fields = fields(contains(fields,vars));
for rm = 1:numel(rm_fields)
  data_avg = rmfield(data_avg,rm_fields{rm});
  if isfield(data_proc,rm_fields{rm})
    data_proc = rmfield(data_proc,rm_fields{rm});
  end
  if isfield(data_proc,[rm_fields{rm} '_rmse'])
    data_proc = rmfield(data_proc,[rm_fields{rm} '_rmse']);
  end
end


%% 5 | Write Level 2 file
keyboard
SUNA_writefile_level2(cfg,data_proc,meta_proc);

%% 5 | Write Level 3 files
keyboard
SUNA_writefile_level3(cfg,data_proc,meta_proc);

%% FUNCTION PLOT_NO3_AND_REFERENCE_DATA
  function plot_NO3_and_reference_data
    makefig; ax=subplot(5,1,1:3); hold(ax,'on'); grid(ax,'on'); ax.YDir = 'normal';  ax.Color = [0.75 0.75 .75]; ax.Box = 'on';
    ax.Title.String = ['SUNA ' cfg.mooring.SN.SUNA ' Reference spectra'];
    ax.YLabel.String = 'NO_3 [\muM]';
    axis(ax,'tight');
    ax2=subplot(5,1,4:5); hold(ax2,'on');grid(ax2,'on');ax2.YDir = 'normal'; ax2.Color = [0.75 0.75 .75]; ax2.Box = 'on';
    ax2.YLabel.String = 'Reference spectra';
    ax2.XLabel.String = 'Wavelength [nm]';
    axis(ax2,'tight');
    ax2.XLim = [190 400];
    % Plot raw and TCSS corrected data with original calibration file
    calstr1 =  ['Cal1: ' calfile_offsets.calfiles{1} ' | ' calfile_offsets.caldates{1} ' | ' num2str(calfile_offsets.lamphour(1),'%.0f') 'lamp hours'];
    refs = struct();
    refs.ref_choices = {calstr1};
    refs.varname     = {'NO3_uM_TCSS_smo'};
    refs.lamphours   = calfile_offsets.lamphour(1);
    refs.lampsec   = calfile_offsets.lampsec(1);
    plot(ax,data_avg.datenum,data_avg.NO3_uM_smo,     'k-','LineWidth',2,'DisplayName',['NO_{3(raw,cal1,smoothed) }  ' calstr1]);
    plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_smo,'g-','LineWidth',2,'DisplayName',['NO_{3(TCSS,cal1,smoothed) } ' calstr1]);
    % Plot TCSS corrected data with calibration file from next deployment
    if isfield(data_avg,'NO3_uM_TCSS2_smo')
      calstr2 =  ['Cal2: ' calfile_offsets.calfiles{2} ' | ' calfile_offsets.caldates{2} ' | ' num2str(calfile_offsets.lamphour(2),'%.0f') 'lamp hours'];
      plot(ax2, calfile_offsets.cal1.Wavelength,calfile_offsets.cal1.Reference,         'g-','LineWidth',2,'DisplayName',calstr1);
      refs.ref_choices = [refs.ref_choices; calstr2];
      refs.varname     = [refs.varname; 'NO3_uM_TCSS2_smo'];
      refs.lamphours   = [refs.lamphours; calfile_offsets.lamphour(2)];
      refs.lampsec     = [refs.lampsec; calfile_offsets.lampsec(2)];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS2_smo,                      'b-','LineWidth',2,'DisplayName',['NO_{3(TCSS,cal2,smoothed) } ' calstr2]);
      plot(ax2,calfile_offsets.cal1.Wavelength,calfile_offsets.cal2.Reference, 'b-','LineWidth',2,'DisplayName',calstr2);
    end
    % Plot TCSS corrected data with pre-deployment reference spectra
    if isfield(data_avg,'NO3_uM_TCSS_pre_smo')
      refs.ref_choices = [refs.ref_choices; [cal_pre.csvfile ' | ' num2str(cal_pre.lamphrs) ' lamp hours']];
      refs.varname     = [refs.varname; 'NO3_uM_TCSS_pre_smo'];
      refs.lamphours   = [refs.lamphours; cal_pre.lamphrs];
      refs.lampsec     = [refs.lampsec;  cal_pre.lampsec];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_pre_smo,      'r-','LineWidth',2,'DisplayName',['NO_{3(TCSS,pre,smoothed) } pre-deployment: ' strrep(cal_pre.csvfile,'_','\_') ' | ' num2str(cal_pre.lamphrs) ' lamp hours' ]);
      plot(ax2,calfile_offsets.cal1.Wavelength,cal_pre.Reference, 'r-','LineWidth',2,'DisplayName',[strrep(cal_pre.csvfile,'_','\_') ' | ' num2str(cal_pre.lamphrs) ' lamp hours']);
    end
    % Plot TCSS corrected data with post-recovery reference spectra
    if isfield(data_avg,'NO3_uM_TCSS_post_smo')
      refs.ref_choices = [refs.ref_choices; [cal_post.csvfile ' | ' num2str(cal_post.lamphrs) ' lamp hours']];
      refs.varname     = [refs.varname; 'NO3_uM_TCSS_post_smo'];
      refs.lampsec     = [refs.lampsec; cal_post.lampsec];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_post_smo,      'm-','LineWidth',2,'DisplayName',['NO_{3(TCSS,post,smoothed) } post-recovery: ' strrep(cal_post.csvfile,'_','\_') ' | ' num2str(cal_post.lamphrs) ' lamp hours' ]);
      plot(ax2,calfile_offsets.cal1.Wavelength,cal_post.Reference, 'm-','LineWidth',2,'DisplayName',[strrep(cal_post.csvfile,'_','\_') ' | ' num2str(cal_post.lamphrs) ' lamp hours']);
    end
    datetick(ax,'x');
    hl = legend(ax,'show');   hl.FontSize = 12;
    hl.Units = 'normalized'; hl.Position(1:2) = [0.002 0.85];
    hl2 = legend(ax2,'show'); hl2.FontSize = 12;
    fprintf('\n');
   
  end %% FUNCTION PLOT_NO3_AND_REFERENCE_DATA

%% FUNCTION PLOT_NO3
  function plot_NO3
    makefig; ax=gca; hold(ax,'on'); grid(ax,'on'); ax.YDir = 'normal';  ax.Color = [0.75 0.75 .75]; ax.Box = 'on';
    ax.Title.String = ['SUNA ' cfg.mooring.SN.SUNA ' Reference spectra'];
    ax.YLabel.String = 'NO_3 [\muM]';
    axis(ax,'tight');
    % Plot raw and TCSS corrected data with original calibration file
    calstr1 =  ['Cal1: ' calfile_offsets.calfiles{1} ' | ' calfile_offsets.caldates{1} ' | ' num2str(calfile_offsets.lamphour(1),'%.0f') 'lamp hours'];
    plot(ax,data_avg.datenum,data_avg.NO3_uM_smo,     'k-','LineWidth',2,'DisplayName',['NO_{3(raw,cal1,smoothed) } ' calstr1 ]);
    plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_smo,'g-','LineWidth',2,'DisplayName',['NO_{3(TCSS,cal1,smoothed) } ' calstr1]);
    % Plot TCSS corrected data with calibration file from next deployment
    if isfield(data_avg,'NO3_uM_TCSS2_smo')
      calstr2 =  ['Cal2: ' calfile_offsets.calfiles{2} ' | ' calfile_offsets.caldates{2} ' | ' num2str(calfile_offsets.lamphour(2),'%.0f') 'lamp hours'];
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS2_smo,    'b-','LineWidth',2,'DisplayName',['NO_{3(TCSS,cal2,smoothed) } ' calstr2]);
    end
    % Plot TCSS corrected data with pre-deployment reference spectra
    if isfield(data_avg,'NO3_uM_TCSS_pre_smo')
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_pre_smo, 'r-','LineWidth',2,'DisplayName',['NO_{3(TCSS,pre,smoothed) } pre-deployment: ' strrep(cal_pre.csvfile,'_','\_') ' | ' num2str(cal_pre.lamphrs) ' lamp hours' ]);
    end
    % Plot TCSS corrected data with post-recovery reference spectra
    if isfield(data_avg,'NO3_uM_TCSS_post_smo')
      plot(ax,data_avg.datenum,data_avg.NO3_uM_TCSS_post_smo,'m-','LineWidth',2,'DisplayName',['NO_{3(TCSS,post,smoothed) } post-recovery: ' strrep(cal_post.csvfile,'_','\_') ' | ' num2str(cal_post.lamphrs) ' lamp hours' ]);
    end
    datetick(ax,'x');
    hl = legend(ax,'show');   hl.FontSize = 12;
    hl.Units = 'normalized'; hl.Position(1:2) = [0.002 0.85];
    fprintf('\n');
  end %% FUNCTION PLOT_NO3


end %% MAIN FUNCTION