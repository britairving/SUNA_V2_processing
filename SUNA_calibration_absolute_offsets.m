function [cfg,data_proc,meta_proc,offsets] = SUNA_calibration_absolute_offsets(data_avg,data_proc,meta_proc,cfg)
%FUNCTION SUNA_calibration_absolute_offsets
%
%  Syntax:
%    [data_proc,meta_proc,cfg] = SUNA_calibration_absolute_offsets(data_proc,meta_proc,cfg)
%
%  Description: Calculate relative offset by comparing TCSS data calculated
%    with different calibration files. Offsets are relative because
%    calibration could still be wrong - so calculated values are not
%    absolute.
%
%  Syntax:
%
%  Description:
%
%  Examples:
%
%  References:
%     Pellerin, B. A., Bergamaschi, B. A., Downing, B. D., Saraceno, J. F.,
%     Garrett, J. D., and Olsen, L. D. (2013). Optical techniques for the
%     determination of nitrate in environmental waters: Guidelines for
%     instrument selection, operation, deployment, maintenance, quality
%     assurance, and data reporting. U.S. Geological Survey Techniques and
%     Methods 1-D5, 37. doi: 10.3133/t m1D5
%
%
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Set up script
if isfield(cfg,'save_figures')
  save_fig = cfg.save_figures;
else
  save_fig = false;
end

lampsec1 = cfg.SUNA.cal1.lampsec;
lamphrs1 = round(cfg.SUNA.cal1.lampsec./60./60,4);
lampstr1 = [num2str(lamphrs1,'%.0f') ' lamp hours'];
caldate1 = cfg.SUNA.cal1.caldate;
calstr1  = cfg.SUNA.cal1.calfile;
if isfield(cfg.SUNA,'cal2')
  lampsec2 = cfg.SUNA.cal2.lampsec;
  lamphrs2 = round(cfg.SUNA.cal2.lampsec./60./60,4);
  lampstr2 = [num2str(lamphrs2,'%.0f') ' lamp hours'];
  caldate2 = cfg.SUNA.cal2.caldate;
  calstr2  = cfg.SUNA.cal2.calfile;
else
  fprintf(' Only OLD calibration files exist..\n')
  offsets = struct();
  offsets.cal1     = cfg.SUNA.cal1.caldata;
  offsets.calfiles = {calstr1};       % Cal file
  offsets.caldates = {caldate1 };         % Date of calibration
  offsets.lamphour = [round(lamphrs1,3)]; % Lamp hours (near) time of calibration
  offsets.lampsec  = [lampsec1];          % Lamp seconds (near) time of calibration
  return
end
try
  %% 1 | Load caldata
  try
    %fprintf('Calculating Nitrate TCSS with second calibration file: %s\n',calfile2)
    %[cal2, ~] = SUNA_read_calfile(calfile2);
    cal2 = cfg.SUNA.cal2.caldata;
  catch
    fprintf('cfg.SUNA.cal2.caldata does not exist... check here\n');
    keyboard
  end
  
  %% 2 | Calculate TCSS data with second calibration file
  if cfg.SUNA.TCSS_ignore_flag  % do not calculate TCSS where data has been flagged
    [NTR,rmse] = ISUS_REPROCESSOR_v2_bi(data_proc,meta_proc,cal2);% pass meta_proc if want to ignore bad data
  else % calculate TCSS for every measurement regardless of flag
    [NTR,rmse] = ISUS_REPROCESSOR_v2_bi(data_proc,[],cal2);       % do NOT pass meta_proc through if want to include bad data
  end
  
  data_proc.NO3_uM_TCSS2 = NTR;
  data_proc.NO3_uM_TCSS2_rmse = rmse;
  % update metadata
  meta_proc.NO3_uM_TCSS2         = meta_proc.NO3_uM;
  meta_proc.NO3_uM_TCSS2.sources = {'NO3_uM' 'dark' 'salinity' 'temperature' 'spectrum_channels'};
  meta_proc.NO3_uM_TCSS2.method  = 'ISUS_REPROCESSOR_v2_bi.m';
  meta_proc.NO3_uM_TCSS2.calfile = cfg.SUNA.cal2.calfile;
  meta_proc.NO3_uM_TCSS2.caldata = cfg.SUNA.cal2.caldata;
  meta_proc.NO3_uM_TCSS2.meaning = 'Sakamoto et al. (2009) Temperature compensated salinity subtracted (TCSS)';
  % meta_proc.NO3_uM_TCSS2_rmse          = struct();
  % meta_proc.NO3_uM_TCSS2_rmse.sources = {'NO3_uM' 'dark' 'salinity' 'temperature' 'spectrum_channels'};
  % meta_proc.NO3_uM_TCSS2_rmse.method  = 'ISUS_REPROCESSOR_v2_bi.m';
  % meta_proc.NO3_uM_TCSS2_rmse.calfile = calfile2;
  % meta_proc.NO3_uM_TCSS2_rmse.caldata = cal2;
  % meta_proc.NO3_uM_TCSS2_rmse.meaning = 'Sakamoto et al. (2009) TCSS RMS error of fit';
  
  %% 3 | Plot Dark current-corrected intensities of DIW/MQ from calibrations
  % All intensities were corrected for dark current using the average
  % intensity of the first 100 pixels (~80 nm wavelength range) measured with
  % the lamp’s shutter closed. The average of only the first 100 pixels is
  % used to eliminate any possible ambient light effects at longer UV
  % wavelengths. [ Sakomoto et al., 2009]
  cal1  = meta_proc.NO3_uM_TCSS.caldata;
  
  %% 4 | Average segment data
  % initialize new variables
  vars = {'NO3_uM_TCSS2' 'NO3_uM_TCSS2_rmse'};
  % Initialize variables
  for nvar = 1:numel(vars)
    vname = vars{nvar};
    data_avg.(vname) = nan(data_avg.number_of_bursts,1);
  end
  % Loop through segments and pull out average counts and nitrate
  for nseg = 1:data_avg.number_of_bursts
    rng = find(data_proc.burst_index == nseg & data_proc.flag <= meta_proc.flag.not_evaluated);
    for nvar = 1:numel(vars)
      vname = vars{nvar};
      data_avg.(vname)(nseg) = mean(data_proc.(vname)(rng));
    end
  end
  
  % smooth
  hr_window = 35;    % 35-hour filter window for movmedian filter
  % Instead.. try movmdedian with 35 hour filter
  filter_vars = {'NO3_uM_TCSS2'}; %{'NO3_uM' 'NO3_uM_TCSS' 'NO3_uM_TCSS2'}
  idx = isfinite(data_avg.datetime);
  for nv = 1:numel(filter_vars)
    smo_var  = [filter_vars{nv} '_smo'];
    % Initialize
    data_avg.(smo_var)      = nan(size(data_avg.(filter_vars{nv})));
    data_avg.(smo_var)(idx) = smoothdata(data_avg.(filter_vars{nv})(idx),'movmedian',hours(hr_window),'SamplePoints',data_avg.datetime(idx));
  end
  
  lamp_hours = data_avg.cum_lampon_sec./60./60;
  makefig; ax=subplot(4,1,1:3); hold(ax,'on'); grid(ax,'on');ax.YDir = 'normal'; ax.Color = [0.75 0.75 .75]; ax.Title.String = [strrep(cfg.project,'_',' ') ' | Lamp degradation'];
  plot(ax,lamp_hours,data_avg.NO3_uM_TCSS_smo,'ko', 'MarkerSize',7,'LineWidth',1,'MarkerFaceColor','g','DisplayName',['NO_{3(TCSS, smoothed)} Cal: ' calstr1 ' on ' caldate1 ' | ' lampstr1 ]);
  plot(ax,lamp_hours,data_avg.NO3_uM_TCSS2_smo,'kd','MarkerSize',7,'LineWidth',1,'MarkerFaceColor','b','DisplayName',['NO_{3(TCSS, smoothed)} Cal: ' calstr2 ' on ' caldate2 ' | ' lampstr2 ]);
  hl = legend(ax,'show'); hl.FontSize = 14; hl.Location = 'southeast';
  ax.YLabel.String = 'NO_3 [\muM]';
  ax2=subplot(4,1,4); hold(ax2,'on'); grid(ax2,'on');ax2.YDir = 'normal'; ax2.Color = [0.75 0.75 .75];
  plot(ax2,lamp_hours,data_avg.NO3_uM_TCSS_smo - data_avg.NO3_uM_TCSS2_smo,'k*', 'MarkerSize',3,'LineWidth',2,'DisplayName','Cal1 - Cal2');
  hl = legend(ax2,'show'); hl.FontSize = 14;hl.Location = 'southeast';
  ax2.XLabel.String = 'Lamp hours';
  ax2.YLabel.String = 'NO_3 Difference [\muM]';
  
  ending_offset = nanmean(data_avg.NO3_uM_TCSS_smo(end-100:end) - data_avg.NO3_uM_TCSS2_smo(end-100:end));
  cal_offset = interp1([lamphrs1 lamphrs2],[0 ending_offset],lamp_hours);
  test = data_avg.NO3_uM_TCSS_smo - cal_offset;
  plot(ax,lamp_hours,test,'ks','MarkerSize',4,'LineWidth',0.2,'MarkerFaceColor','y','DisplayName','NO_{3(TCSS, smoothed)} corrected with cal2 offset');
  plot(ax2,lamp_hours,cal_offset,'y*', 'MarkerSize',3,'LineWidth',2,'DisplayName','Interpolated offsets calculated with cal2');
  ax2.XLim = ax.XLim;
  
  if save_fig
    savename = fullfile(cfg.datadir,[cfg.project '_calibration_comparison_lampdegradation1']);
    standard_printfig_highrespng(savename);
  end
  
  %% 5 | Update offsets structure to return to SUNA_correct_data.m script
  offsets = struct();
  offsets.cal1     = cal1; 
  offsets.cal2     = cal2;
  offsets.calfiles = {calstr1  calstr2};                    % Cal file
  offsets.caldates = {caldate1 caldate2};                   % Date of calibration
  offsets.lamphour = [round(lamphrs1,3) round(lamphrs2,3)]; % Lamp hours (near) time of calibration
  offsets.lampsec  = [lampsec1 lampsec2];                   % Lamp seconds (near) time of calibration
  offsets.aboffset = [0 ending_offset];                     % Absolute offset calculated with different calibration files
  
  %% 6 | Plot difference reference
  makefig; ax=subplot(3,1,1:2); hold(ax,'on'); grid(ax,'on');ax.YDir = 'normal'; ax.Color = [0.75 0.75 .75]; ax.Title.String = ['SUNA ' cfg.mooring.SN.SUNA ' Reference spectra'];
  plot(cal1.Wavelength,cal1.Reference,'ko', 'MarkerSize',7,'LineWidth',1,'MarkerFaceColor','g','DisplayName',['Cal: ' calstr1 ' on ' caldate1 ' | ' lampstr1 ]);
  plot(cal2.Wavelength,cal2.Reference,'kd','MarkerSize',7,'LineWidth',1,'MarkerFaceColor','b','DisplayName',['Cal: ' calstr2 ' on ' caldate2 ' | ' lampstr2 ]);
  hl = legend(ax,'show'); hl.FontSize = 14;
  ax.YLabel.String = 'Reference Intensity';
  
  ax.XLim = [190 400];
  ax2=subplot(3,1,3); hold(ax2,'on'); grid(ax2,'on');ax2.YDir = 'normal'; ax2.Color = [0.75 0.75 .75];
  % Percent difference
  pdiff = (cal1.Reference - cal2.Reference)./cal1.Reference * 100.0;
  idx = cal1.Wavelength > 190;
  plot(ax2,cal1.Wavelength(idx),pdiff(idx),'k*', 'MarkerSize',7,'LineWidth',3,'DisplayName','% difference');
  ax2.XLim = ax.XLim;
  ax2.XLabel.String = 'Wavelength [nm]';
  ax2.YLabel.String = 'difference (%)';
  
  if save_fig
    savename = fullfile(cfg.datadir,[cfg.project '_calibration_comparison_lampdegradation2']);
    standard_printfig_highrespng(savename);
  end
  
  %% 7 | save the data
  calfile_offsets = offsets;
  save(cfg.path.calfile_step1,'data_proc','meta_proc','calfile_offsets');
catch
  fprintf('Some error in SUNA_calibration_absolute_offsets.m -- investigate further\n')
  fprintf('Enter keyboard mode... enter "dbcont" to continue\n');
  keyboard
end
end %% MAIN FUNCTION