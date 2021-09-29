function [cfg,data_pre,meta_pre] = SUNA_preprocessing(cfg,data_raw,meta_raw)
%FUNCTION SUNA_preprocessing
%
%  Syntax:
%    [cfg, t] = SUNA_preprocessing(cfg,data_raw)
%
%  Description:
%    SeapHOx_process_data
%
%  Syntax:
%    [DATA_PRE, META_PRE] = PREPROCESSGLIDERDATA(DATA_RAW)
%    [DATA_PRE, META_PRE] = PREPROCESSGLIDERDATA(DATA_RAW, OPTIONS)
%    [DATA_PRE, META_PRE] = PREPROCESSGLIDERDATA(DATA_RAW, OPT1, VAL1, ...)
%
%  Description:
%    [DATA_PRE, META_PRE] = PREPROCESSGLIDERDATA(DATA_RAW, ...)
%    preprocesses seapHOx deployment data according to given options,
%    performing the following actions:
%      - Selection of horizontal position sensors:
%        Latitude and longitude sequences are selected.
%        A position flag is also selected, if any, and bad values are masked
%        as missing. Optionally, a unit conversion may be applied, if needed.
%        Latitude and longitude are mandatory, processing aborts if missing.
%      - Selection of optional reference sensors:
%        Navigation depth, pitch, roll, and heading sequences are selected,
%        if any. Unit conversions may be applied, if needed.
%      - Selection of water velocity sensors or estimates.
%        Mean segment water eastward and northward velocity sequences are
%        selected, if any. Unit conversions may be applied, if needed.
%      - Selection of commanded trajectory sensors:
%        Commanded waypoint longitude and latitude sequences are selected,
%        if any. Unit conversions may be applied, if needed.
%      - Selection of CTD sensor:
%        Conductivity, temperature and pressure sequences are selected, if any.
%        Optionally the CTD timestamp sequence is selected, too.
%        Unit conversion may be applied to pressure readings, if needed; and
%        factory calibrations may be applied to temperature and conductivity.
%      - Selection of fluorescence and/or scattering sensor:
%        Chlorophyll, turbidity, CDOM and/or scattering sequences are selected,
%        if any. Optionally the sensor timestamp sequence is selected, too.
%        Manufacturer calibrations may be applied, if needed.
%      - Selection of oxygen sensors:
%        Oxygen concentration and saturation sequences are selected, if any.
%        Optionally oxygen sensor timestamp and temperature sequences are
%        selected. Unit conversion and manufacturer calibrations may be applied.
%      - Selection of other sensors of interest:
%        Sequences from extra sensors configured in options are selected.
%        Unit conversions and manufacturer calibrations may be applied,
%        if needed.
%
%
%  Examples:
%    data_pre = preprocessGliderData(data_raw, options)
%
%  References:
%    Martz, T., K. McLaughlin, S.B. Weisberg. 2015. Best Practices for
%    autonomous measurement of seawater pH with the Honeywell Durafet pH
%    sensor. California Current Acidification Network (C-CAN).
%
%    Bresnahan et al., 2014. P.J. Bresnahan, T.R. Martz, Y. Takeshita, K.S.
%    Johnson, M. LaShomb Best practices for autonomous measurement of
%    seawater pH with the Honeywell Durafet Methods Oceanogr., 9 (2014),
%    pp. 44-60
%
%    National Academies of Sciences, Engineering, and Medicine 2017.
%    Valuing Climate Damages: Updating Estimation of the Social Cost of
%    Carbon Dioxide. Washington, DC: The National Academies Press.
%    https://doi.org/10.17226/24651.
%
%    Takahashi, T., S.C. Sutherland, D.W. Chipman, J. G. Goddard, C. Ho, T.
%    Newberger, C. Sweeney, and D. R. Munro. 2014. Climatological
%    distributions of pH, pCO2, total CO2, alkalinity, and CaCO3 saturation
%    in the global surface ocean, and temporal changes at selected
%    locations. Mar. Chem. 164: 95–125. doi:10.1016/j.marchem.2014. 06.004
%
%    Garcia, H. E., R. A. Locarnini, T. P. Boyer, J. I. Antonov, M. M.
%    Zweng, O. K. Baranova, and D. R. Johnson, 2010. World Ocean Atlas
%    2009, Volume 4: Nutrients (phosphate, nitrate, silicate). S. Levitus,
%    Ed. NOAA Atlas NESDIS 71, U.S. Government Printing Office, Washington,
%    D.C., 398 pp.
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Script set up
dbstop if error

%% 1 | Initialize data_pre and meta_pre structures
% meta_pre will have all the same fieldnames as data_pre and store
% information about variable such as units, conversions, etc.
data_pre = struct();
meta_pre = struct();

data_pre.project  = cellstr(repmat(cfg.project,size(data_raw.yyyydoy)));
data_pre.datetime = datetime(string(data_raw.yyyydoy),'InputFormat','yyyyDDD') + hours(data_raw.hr);
data_pre.datenum  = datenum(data_pre.datetime);
meta_pre.project        = data_pre.project{1};
meta_pre.datetime.units = 'DateTime';
meta_pre.datenum.units  = 'Days since January 0, 0000';

data_pre.latitude  = repmat(cfg.mooring.latitude,size(data_pre.datenum));
data_pre.longitude = repmat(cfg.mooring.longitude,size(data_pre.datenum));
meta_pre.latitude.units = 'degrees north';
meta_pre.longitude.units = 'degrees east';
try
  %% 1 | Pull out relevant variables from raw data
  data_pre.cum_lampon_sec = data_raw.cum_lampon_sec;
  meta_pre.cum_lampon_sec = meta_raw.cum_lampon_sec;
  
  data_pre.fit_RMSE = data_raw.fit_RMSE;
  meta_pre.fit_RMSE = meta_raw.fit_RMSE;
  meta_pre.fit_RMSE.meaning = 'The root mean square error parameter from the SUNA V2 can be used to make an estimate of how well the nitrate spectral fit is. This should usually be less than 1E-3. If it is higher, there is spectral shape (likely due to CDOM) that adversely impacts the nitrate estimate.';
  
  data_pre.abs_254nm = data_raw.abs_254nm;
  meta_pre.abs_254nm = meta_raw.abs_254nm;
  meta_pre.abs_254nm.long_name = 'Absorbance at 254nm';
  meta_pre.abs_254nm.meaning   = 'Can be used to make an estimate of the impact of CDOM. If absorption is high (> 1.3 AU), the SUNA will not be able to collect adequate light to make a measurement.';
  
  data_pre.abs_350nm = data_raw.abs_350nm;
  meta_pre.abs_350nm = meta_raw.abs_350nm;
  meta_pre.abs_350nm.long_name = 'Absorbance at 350nm';
  meta_pre.abs_350nm.meaning   = 'Can be used to make an estimate of the impact of CDOM. If absorption is high (> 1.3 AU), the SUNA will not be able to collect adequate light to make a measurement.';
  
  data_pre.bromide_mgL = data_raw.bromide_mgL;
  meta_pre.bromide_mgL = meta_raw.bromide_mgL;
  meta_pre.bromide_mgL.long_name = 'Bromide trace [mg/L]';
  
  data_pre.dark      = data_raw.dark;
  meta_pre.dark.units = 'Dark value used for fit';
  
  data_pre.NO3_uM      = data_raw.NO3_uM;
  meta_pre.NO3_uM.units = 'umol/L';
  meta_pre.NO3_uM.long_name     = 'Mole Concentration Of Nitrate In Sea Water';
  meta_pre.NO3_uM.standard_name = 'mole_concentration_of_nitrate_in_sea_water';
  meta_pre.NO3_uM.ioos_category = 'Dissolved Nutrients';
  meta_pre.NO3_uM.comment       = '1 µM nitrogen = 1 umol/L nitrogen = 0.014007 mg/L nitrogen';
  
  data_pre.NO3_mgNL      = data_raw.NO3_mgNL;
  meta_pre.NO3_mgNL.units  = 'mgN/L';
  meta_pre.NO3_mgNL.method = 'converted from NO3_uM where 1 µM nitrogen = 0.014007 mg/L nitrogen';

  % spectrum channels
  fields = fieldnames(data_raw);
  spec_fields = contains(fields,'spec') & ~contains(fields,{'spec_temp' 'spec_avg'});
  data_pre.spectrum_channels = table2array(data_raw(:,spec_fields));
  meta_pre.spectrum_channels.units = 'raw spectrum, channels 1-256';
catch
  fprintf(' problem here assigning data_pre variables from data_raw table..\n')
  keyboard
end

%% 1 | Read in Level 1 CTD data
% Overwritten with Level2 data below.. just read in level1 because it
% contains pre/post deployment data which are thrown out in level 2 and in
% case it is higher resolution
if isfield(cfg.SUNA.ctd,'level1') && ~isempty(cfg.SUNA.ctd.level1)
  cfg.SUNA.ctd.file  = cfg.SUNA.ctd.level1;
  cfg.SUNA.ctd.level = 'Level 1';
  try
    ctd = cnv2mat_mooring(cfg.SUNA.ctd.file);
  catch
    keyboard
    ctd_opt = detectImportOptions(cfg.SUNA.ctd.file,'FileType','text','CommentStyle', {'*' '#'});
    ctd_opt.CommentStyle = {'*' '#'};
  end
end

%% 2 | Calculate time lag 
% This just stops if user decides there is a time lag, have not actually
% implemented the time lag correction because have not encountered yet. 
makefig; ax = gca; ax.YDir = 'normal'; hold(ax,'on');grid(ax,'on');
gd = data_pre.dark ~= 0;
plot(ax,data_pre.datenum(gd),data_pre.NO3_uM(gd),'b*','LineWidth',2,'DisplayName','NO_3 uM')
plot(ax,ctd.datenum,ctd.salinity,'r*','LineWidth',2,'DisplayName','Salinity')
plot(ax,ctd.datenum,ctd.pressure,'ks','LineWidth',2,'DisplayName','Pressure')
datetick(ax,'x','mm/yyyy','Keepticks','KeepLimits');
ax.XTickLabelRotation = 45;
ax.XAxis.FontSize = 14;
ax.YLim = [-100 100];
try
  pause(0.1);
  hl = legend(ax,'show');
  hl = legend(ax,'show');
catch
  text(ax,0.9,0.98,'NO_3 uM','FontSize',16,'Color','b','Units','Normalized','FontWeight','bold')
  text(ax,0.9,0.94,'Salinity','FontSize',16,'Color','r','Units','Normalized','FontWeight','bold')
  text(ax,0.9,0.90,'Pressure','FontSize',16,'Color','k','Units','Normalized','FontWeight','bold')
end
text(ax,0.01,0.98,'Zoom in to examine if there is a timelag between the CTD data and the SUNA data','FontSize',16,'Color','r','Units','Normalized','FontWeight','bold')
text(ax,0.01,0.95,'For example... if there is a clear lag between when values change at deployment/recovery','FontSize',16,'Color','r','Units','Normalized','FontWeight','bold')

tlag_chc = input('Is there a time lag? <1/0> ');
if isempty(tlag_chc); tlag_chc = 0; end % default to NO
if tlag_chc == 1
  fprintf('Need to write up method to correct time lag here!!\n')
  keyboard
end

%% 1 | Read in Level 2 CTD data, if available
if isfield(cfg.SUNA.ctd,'level2') && ~isempty(cfg.SUNA.ctd.level2)
  cfg.SUNA.ctd.file  = cfg.SUNA.ctd.level2;
  cfg.SUNA.ctd.level = 'Level 2';
  ctd = read_level2_ctdfile(cfg.SUNA.ctd.file);
end

%% 1 | CTD Data
% Add CTD data
if all(isnan(data_raw.ctd_temp)) && exist(cfg.SUNA.ctd.file,'file')
  % Pull out ctd filename without folder to store
  [~,ctd_filename,ctd_fileext] = fileparts(cfg.SUNA.ctd.file);
  ctd_filename = [ctd_filename ctd_fileext];
  
  % Interpolate time
  if isfield(cfg.SUNA.ctd,'SUNA_timelag')
    fprintf('add time lag correciton here\n')
    keyboard
  end
  % Pressure
  data_pre.pressure         = interp1(ctd.datenum,ctd.pressure,data_pre.datenum);
  meta_pre.pressure.sources = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
  meta_pre.pressure.units    = 'dbar';
  meta_pre.pressure.ctdfile = ctd_filename;
  meta_pre.pressure.standard_name = 'sea_water_pressure';
  meta_pre.pressure.instrument    = 'instrument_ctd';
  if isfield(ctd,'pressure_flag')
    data_pre.pressure_flag         = interp1(ctd.datenum,ctd.pressure_flag,data_pre.datenum);
    meta_pre.pressure_flag.sources = ctd_filename;
  end
  % Temperature
  data_pre.temperature         = interp1(ctd.datenum,ctd.temp,data_pre.datenum);
  meta_pre.temperature.sources = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
  meta_pre.temperature.units    = 'degC';
  meta_pre.temperature.ctdfile = ctd_filename;
  meta_pre.temperature.standard_name = 'sea_water_temperature';
  meta_pre.temperature.instrument    = 'instrument_ctd';
  if isfield(ctd,'temp_flag')
    data_pre.temp_flag         = interp1(ctd.datenum,ctd.temp_flag,data_pre.datenum);
    meta_pre.temp_flag.sources = ctd_filename;
  end
  % Conductivity
  if isfield(ctd,'cond')
    data_pre.conductivity         = interp1(ctd.datenum,ctd.cond,data_pre.datenum);
    meta_pre.conductivity.sources = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    meta_pre.conductivity.units   = ctd.units.cond;
    meta_pre.conductivity.ctdfile = ctd_filename;
    meta_pre.conductivity.standard_name = 'sea_water_electrical_conductivity';
    meta_pre.conductivity.instrument    = 'instrument_ctd';
  end
  if isfield(ctd,'cond_flag')
    data_pre.conductivity_flag         = interp1(ctd.datenum,ctd.cond_flag,data_pre.datenum);
    meta_pre.conductivity_flag.sources = ctd_filename;
  end
  % Salinity
  if isfield(ctd,'salinity')
    data_pre.salinity         = interp1(ctd.datenum,ctd.salinity,data_pre.datenum);
    meta_pre.salinity.sources = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    meta_pre.salinity.ctdfile = ctd_filename;
  else
    fprintf('This does not calculate salinity.. just sets it up. Need to finish here\n')
    keyboard
    data_pre.salinity         = interp1(ctd.datenum,gsw_SP_from_C(ctd.cond*10,ctd.temp,ctd.pressure),data_pre.datenum);
    meta_pre.salinity.sources = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    meta_pre.salinity.method  = 'gsw_SP_from_C';
  end
  if isfield(ctd,'salinity_flag')
    data_pre.salinity_flag         = interp1(ctd.datenum,ctd.salinity_flag,data_pre.datenum);
    meta_pre.salinity_flag.sources = ctd_filename;
  end
  if isfield(ctd.units,'salinity')
    meta_pre.salinity.units  = ctd.units.salinity;
  else
    meta_pre.salinity.units  = 'unitless';
  end  
  meta_pre.salinity.standard_name    = 'sea_water_salinity';
  meta_pre.pressure.standard_name    = 'sea_water_pressure';
  meta_pre.temperature.standard_name = 'sea_water_temperature';
  % Calculate density (TEOS 10)
  if isfield(data_pre,'salinity') && isfield(data_pre,'temperature') && isfield(data_pre,'pressure')
    % Use gsw_rho_irving because automatically calculates necessary
    % variables for gsw_rho
    d.p = data_pre.pressure;
    d.t = data_pre.temperature;
    d.Sp = data_pre.salinity;
    d.lat = data_pre.latitude;
    d.lon = data_pre.longitude;
    dens = gsw_rho_irving(d);
    % store in-situ density of seawater
    data_pre.density = dens.rho;
    meta_pre.density.sources = {'pressure' 'salinity' 'temperature' 'latitude' 'longitude'};
    meta_pre.density.method  = 'gsw_rho_irving';
    meta_pre.density.units   = 'kg/m3';
    meta_pre.density.standard_name = 'sea_water_density';
    meta_pre.density.instrument    = 'instrument_ctd';
    % store potential density
    data_pre.potential_density = dens.pot_rho;
    meta_pre.potential_density.sources = {'pressure' 'salinity' 'temperature' 'latitude' 'longitude'};
    meta_pre.potential_density.method  = 'gsw_rho_irving';
    meta_pre.potential_density.units   = 'kg/m3';
    meta_pre.potential_density.standard_name = 'sea_water_potential_density';
    meta_pre.potential_density.instrument    = 'instrument_ctd';
  end
  %% ECO sensors
  % CDOM
  if isfield(ctd,'cdom')
    data_pre.cdom          = interp1(ctd.datenum,ctd.cdom,data_pre.datenum);
    meta_pre.cdom.sources  = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    meta_pre.cdom.standard_name = 'concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent_mass_fraction_of_quinine_sulfate_dihydrate';
    if isfield(ctd.units,'cdom')
      meta_pre.cdom.units = ctd.units.cdom;
    else
      meta_pre.cdom.units    = ctd.units{contains(ctd.names,'cdom')};
      if contains( meta_pre.cdom.units,'[mg/m^3]')
        meta_pre.cdom.long_name = meta_pre.cdom.units;
        meta_pre.cdom.units = 'mg/m^3';
      else
        fprintf('check CDOM units and try to simplify\n')
        keyboard
      end
    end
    meta_pre.cdom.ctdfile  = ctd_filename;
  end
  % Fluorometer
  if isfield(ctd,'fluor_mgm3')
    data_pre.fluor_mgm3         = interp1(ctd.datenum,ctd.fluor_mgm3,data_pre.datenum);
    meta_pre.fluor_mgm3.sources = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    if isfield(ctd.units,'fluor_mgm3')
      meta_pre.fluor_mgm3.units = ctd.units.fluor_mgm3;
    else
      meta_pre.fluor_mgm3.units   = ctd.units{contains(ctd.names,'fluor_mgm3')};
      if contains( meta_pre.fluor_mgm3.units,'[mg/m^3]')
        meta_pre.fluor_mgm3.long_name = meta_pre.fluor_mgm3.units;
        meta_pre.fluor_mgm3.units = 'mg/m^3';
      else
        fprintf('check fluor_mgm3 units and try to simplify\n')
        keyboard
      end
    end
    meta_pre.fluor_mgm3.ctdfile = ctd_filename;
  end
  % turbidity
  if isfield(ctd,'turbidity')
    data_pre.turbidity          = interp1(ctd.datenum,ctd.turbidity,data_pre.datenum);
    meta_pre.turbidity.sources  = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    meta_pre.turbidity.standard_name = 'sea_water_turbidity';
    meta_pre.turbidity.instrument    = 'instrument_ctd';
    if isfield(ctd.units,'turbidity')
      meta_pre.turbidity.units = ctd.units.turbidity;
    else
      meta_pre.turbidity.units    = ctd.units{contains(ctd.names,'turbidity')};
      if contains( meta_pre.turbidity.units,'[m^-1/sr]')
        meta_pre.turbidity.long_name = meta_pre.turbidity.units;
        meta_pre.turbidity.units = 'm^-1/sr';
      else
        fprintf('check turbidity units and try to simplify\n')
        keyboard
      end
    end
    meta_pre.turbidity.ctdfile  = ctd_filename;
  end
  % PAR
  if isfield(ctd,'PAR')
    data_pre.PAR          = interp1(ctd.datenum,ctd.PAR,data_pre.datenum);
    meta_pre.PAR.sources  = [cfg.SUNA.ctd.type ' SN:' cfg.SUNA.ctd.SN];
    if isfield(ctd.units,'PAR')
      meta_pre.PAR.units = ctd.units.PAR;
    else
      meta_pre.PAR.units    = ctd.units{contains(ctd.names,'PAR')};
      if contains(meta_pre.PAR.units,'Biospherical/Licor')
        meta_pre.PAR.long_name = meta_pre.PAR.units;
        meta_pre.PAR.units = 'Biospherical/Licor';
      else
        fprintf('check PAR units and try to simplify\n')
        keyboard
      end
    end
    meta_pre.PAR.ctdfile  = ctd_filename;
  end
else  %% ~ [  all(isnan(data_raw.ctd_temp)) && exist(cfg.SUNA.ctd.file,'file') ]
  if all(isnan(data_raw.ctd_temp)) && ~exist(cfg.SUNA.ctd.file,'file')
    fprintf('CTD data not avaiable in SUNA files but external ctd file not found: %s\n',cfg.SUNA.ctd.file)
    keyboard
  else
    fprintf('unrecognized scenario... ctd data is in SUNA file??\n')
    keyboard
  end
end

if ~isfield(data_pre,'depth')
  data_pre.depth = -gsw_z_from_p(data_pre.pressure,data_pre.latitude);
  meta_pre.depth.units = 'm';
  meta_pre.depth.method = 'gsw_z_from_p';
end

%% | Diagnostic plots
%% Spectrum plot through time
% fields = fieldnames(data_raw);
% spec_fields = fields(contains(fields,'spec') & ~contains(fields,{'spec_temp' 'spec_avg'}));
% makefig; ax = gca; hold(ax,'on');
% for ns = 1:numel(spec_fields)
%   spec_field = spec_fields{ns};
%   plot(data_pre.datetime,data_raw.(spec_field))
% end
% ax.Title.String = [strrep(cfg.project,'_',' ') ' SUNA Spectrum'];
% grid(ax,'on');
% ax.YLabel.String = 'Raw Spectrum channels 1-256';
% if savefig
%   standard_printfig_highrespng(save_figures,strrep(ax.Title.String,' ','_'))
% end
% % keyboard


% % %% Calibrations for SUNA
% % % default to no
% % precal  = 0;
% % postcal = 0;
% % if isfield(cfg.SUNA.cal,'precal')
% %   precal = 1;
% % end
% % if isfield(cfg.SUNA.cal,'postcal')
% %   postcal = 1;
% % end
% %
% %   % Case 1: precal and postcal measurements were collected
% % if precal && postcal
% %   if ischar(cfg.SUNA.cal.postcal.files)
% %     cfg.SUNA.cal.postcal.files =  cellstr(cfg.SUNA.cal.postcal.files);
% %   end
% %   if ischar(cfg.SUNA.cal.precal.files)
% %     cfg.SUNA.cal.precal.files =  cellstr(cfg.SUNA.cal.precal.files);
% %   end
% %   % build precal dir
% %   caldir = fullfile(cfg.datadir,'CAL');
% %   % PRECAL or PRE-DEPLOYMENT measurements
% %   for nprecal = 1:numel(cfg.SUNA.cal.precal.files)
% %     calfile = cfg.SUNA.cal.precal.files{nprecal};
% %     %  build full filename with path
% %     if exist(fullfile(caldir,calfile),'file')
% %       calfile = fullfile(caldir,calfile);
% %     elseif exist(fullfile(caldir,'deployment',calfile),'file')
% %       calfile = fullfile(caldir,'deployment',calfile);
% %     else
% %       fprintf('could not find %s...\n',calfile)
% %       keyboard
% %     end
% %     % read in data
% %     [~,data_cal,meta_cal] = SUNA_read_data(cfg,calfile);
% %     if nprecal == 1
% %       data_precal = data_cal;
% %     else
% %       data_precal = [data_precal; data_cal];
% %     end
% %   end
% %   % add Matlab's datetime and datenum
% %   data_precal.datetime = datetime(string(data_precal.yyyydoy),'InputFormat','yyyyDDD') + hours(data_precal.hr);
% %   data_precal.datenum  = datenum(data_precal.datetime);
% %
% %   % POSTCAL or POST-DEPLOYMENT measurements
% %   for npostcal = 1:numel(cfg.SUNA.cal.postcal.files)
% %     calfile = cfg.SUNA.cal.postcal.files{npostcal};
% %     %  build full filename with path
% %     if exist(fullfile(caldir,calfile),'file')
% %       calfile = fullfile(caldir,calfile);
% %     elseif exist(fullfile(caldir,'deployment',calfile),'file')
% %       calfile = fullfile(caldir,'deployment',calfile);
% %     else
% %       fprintf('could not find %s...\n',calfile)
% %       keyboard
% %     end
% %     % read in data
% %     [~,data_cal,meta_cal] = SUNA_read_data(cfg,calfile);
% %     if npostcal == 1
% %       data_postcal = data_cal;
% %     else
% %       data_postcal = [data_postcal; data_cal];
% %     end
% %   end
% %   % add Matlab's datetime and datenum
% %   data_postcal.datetime = datetime(string(data_postcal.yyyydoy),'InputFormat','yyyyDDD') + hours(data_postcal.hr);
% %   data_postcal.datenum  = datenum(data_postcal.datetime);
% %   % Diagnostic plots ------------------------------------------------------
% %   makefig; ax = gca;
% %   plot(data_precal.datetime,data_precal.NO3_uM,'b+','LineWidth',1,'MarkerSize',10,'DisplayName','Pre-Deployment')
% %   hold(ax,'on'); ylabel(ax,'Nitrate [\muM]')
% %   plot(data_postcal.datetime,data_postcal.NO3_uM,'bo','LineWidth',1,'DisplayName','Post-Deployment')
% %   grid(ax,'on'); ymin = ax.YLim(1); ymax = ax.YLim(2);
% %   try % try showing lines with deployment and recovery dates
% %     plot(ax,[cfg.mooring.deploydate cfg.mooring.deploydate],[ymin ymax],'k--','DisplayName','Deployment')
% %     plot(ax,[cfg.mooring.recoverdate cfg.mooring.recoverdate],[ymin ymax],'k:','LineWidth',2,'DisplayName','Recovery')
% %   catch; end
% %   plot(ax,data_pre.datetime,data_raw.NO3_uM,'g.','DisplayName','Raw data')
% %   ax.YLim = [ymin ymax];
% %   legend(ax,'Show','Location','north');
% %   if ~isnan(cfg.SUNA.cal.precal.NO3_uM)
% %     plot(ax,data_precal.datetime(1),cfg.SUNA.cal.precal.NO3_uM,'r+','LineWidth',1,'MarkerSize',10,'DisplayName','Pre-Deployment Nitrate concentration')
% %   end
% %   if ~isnan(cfg.SUNA.cal.postcal.NO3_uM)
% %     plot(ax,data_postcal.datetime(1),cfg.SUNA.cal.postcal.NO3_uM,'ro','LineWidth',1,'DisplayName','Post-Deployment Nitrate concentration')
% %   end
% %   % Diagnostic plots ------------------------------------------------------
% %   keyboard
% % % Case 2: ONLY postcal measurements were collected
% % elseif ~precal && postcal
% %   keyboard
% %
% % else
% %   fprintf('set up this scenario....\n')
% %   keyboard
% % end
% % keyboard
% % %% 6 | Save data_pre
% % fprintf('  saving preprocessed data to %s\n',cfg.SUNA.path.data_pre)
% % save(cfg.SUNA.path.data_pre,'data_raw','meta_raw','data_pre','meta_pre','cfg','-v7.3');

end %% function [cfg, t] = SeapHOx_find_zero(cfg,data)