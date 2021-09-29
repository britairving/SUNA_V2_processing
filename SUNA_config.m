function cfg = SUNA_config(cfg)
% FUNCTION SUNA_CONFIG
%
%  Syntax:
%    cfg = SUNA_config(cfg)
%  
%  Description:
%    Build initial metadata by outlining project details, what variables,
%    calibration files, processing steps, etc etc
%
%  See also:
%    User manual Submersible Ultraviolet Nitrate Analyzer SUNA V2, 2017
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% processing flags
cfg.SUNA.TCSS_ignore_flag = 0; % Ignores flagged data when calculating TCSS corrected data

%% Define matlab data filenames
cfg.SUNA.path.data_raw  = fullfile(cfg.datadir,[cfg.project '_SUNA_raw.mat']);  % LEVEL1: raw
cfg.SUNA.path.data_pre  = fullfile(cfg.datadir,[cfg.project '_SUNA_pre.mat']);  % pre processed - converted
cfg.SUNA.path.data_proc = fullfile(cfg.datadir,[cfg.project '_SUNA_proc.mat']); % LEVEL3: processed data
cfg.SUNA.path.data_cor  = fullfile(cfg.datadir,[cfg.project '_SUNA_cor.mat']);  % LEVEL4: calibration drift and/or offset corrected
cfg.SUNA.path.data_qc   = fullfile(cfg.datadir,[cfg.project '_SUNA_qc.mat']);   % LEVEL4: QC'd data
%% Set instrument nominal depth - used to flag data before/after deployment
if contains(cfg.project,'CEO')
  cfg.SUNA.depth = 33;
else
  error('SUNA nominal depth not set yet... set here\n')
  keyboard
end

%% Set calibration and CTD information
switch cfg.project
  case 'CEO_2015'
    % calibration information and files
    cfg.SUNA.cal.postcal = struct();
    cfg.SUNA.cal.postcal.files  = {'SUNA-0562_2016-08-08_20-26-03_raw.csv' 'SUNA-0562_2016-08-08_20-33-51_raw.csv' 'SUNA-0562_2016-08-08_20-47-03_raw.csv'};
    cfg.SUNA.cal.postcal.ctd_file   = '';  % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.postcal.NO3_uM     = NaN; % Nitrate standard concentration in microM. Set to NaN if unknown - defaults to 30uM
    % ctd information
    cfg.SUNA.ctd.type = 'SBE16';
    cfg.SUNA.ctd.SN   = cfg.mooring.SN.(cfg.SUNA.ctd.type);
    %cfg.SUNA.ctd.file = 'F:\CEO\DataByDeploymentYear\2015\SBE16\process\CEM2_2015_SBE16plus_6670_2016_08_05.cnv';
    % This should be named CEO2_2015_SBE16_6670 ...
    cfg.SUNA.ctd.level1 = fullfile(cfg.yeardir,'SBE16','CEO1_2015_SBE16_6670.asc');         % Full resolution- use to determine time lag
    cfg.SUNA.ctd.level2 = fullfile(cfg.yeardir,'SBE16','CEO2_2015-2016_ctd033m_L2_v2.csv'); % QC'd, averaged
    % Calibration file information - created BEFORE deployment
    cfg.SUNA.cal1.calfile = 'SNA0562D.cal';
    cfg.SUNA.cal1.caldate = '16-Aug-2015';
    cfg.SUNA.cal1.lampsec = 184270;    %cum_lampon_sec when (or very near when) cal file was created
    % Calibration file information - created AFTER recovery
    cfg.SUNA.cal2.calfile = 'SNA0562J.cal';
    cfg.SUNA.cal2.caldate = '03-Aug-2017';
    cfg.SUNA.cal2.lampsec = 309390;    %cum_lampon_sec when (or very near when) cal file was created
    
  case 'CEO_2016'
    % calibration information and files
    % predeployment measurements
    cfg.SUNA.cal.precal = struct();
    cfg.SUNA.cal.precal.files   = {'SUNA-0801_2016-08-02_21-40-00_raw.csv' 'SUNA-0801_2016-08-02_21-55-20_raw.csv' 'SUNA-0801_2016-08-02_22-27-03_raw.csv'};
    cfg.SUNA.cal.precal.ctd_file   = '';    % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.precal.NO3_uM     = NaN;   % Nitrate standard concentration in uM. Set to NaN if unknown - defaults to 30uM
    % post deployment measurements
    cfg.SUNA.cal.postcal = struct();
    cfg.SUNA.cal.postcal.files  = 'D2017229_di_then_nuts_postcal.CSV';
    cfg.SUNA.cal.postcal.ctd_file   = '';    % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.postcal.NO3_uM     = NaN;   % Nitrate standard concentration in uM. Set to NaN if unknown - defaults to 30uM
     % ctd information
    cfg.SUNA.ctd.type  = 'SBE16';
    cfg.SUNA.ctd.SN    = cfg.mooring.SN.(cfg.SUNA.ctd.type);
    cfg.SUNA.ctd.level1 = fullfile(cfg.yeardir,'SBE16','CEO2_2016-2017_SBE16_4604.cnv');   % Full resolution- use to determine time lag
    cfg.SUNA.ctd.level2 = fullfile(cfg.yeardir,'SBE16','CEO2_2016-2017_ctd033m_L2_v2.csv');% QC'd, averaged
    % Calibration file information - created BEFORE deployment
    cfg.SUNA.cal1.calfile = 'SNA0801D.cal';
    cfg.SUNA.cal1.caldate = '02-Aug-2016';
    cfg.SUNA.cal1.lampsec = 190380;    %cum_lampon_sec when (or very near when) cal file was created
    % Calibration file information - created AFTER recovery
    cfg.SUNA.cal2.calfile = 'SNA0801E.cal';
    cfg.SUNA.cal2.caldate = '24-Jul-2018';
    cfg.SUNA.cal2.lampsec = 416500;    %cum_lampon_sec when (or very near when) cal file was created
    
  case 'CEO_2017'
    % calibration information and files
    cfg.SUNA.cal.postcal = struct();
    cfg.SUNA.cal.postcal.files  = 'D2018221_contains_DIandStandard_measurements.CSV';
    cfg.SUNA.cal.postcal.ctd_file   = '';  % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.postcal.NO3_uM     = 30;  % Nitrate standard concentration in uM. Set to NaN if unknown - defaults to 30uM
    % ctd information
    cfg.SUNA.ctd.type = 'SBE16';
    cfg.SUNA.ctd.SN   = cfg.mooring.SN.(cfg.SUNA.ctd.type);
    cfg.SUNA.ctd.level1 = fullfile(cfg.yeardir,'SBE16','CEO2_2017-2018_SBE16_6670.cnv');   % Full resolution- use to determine time lag
    cfg.SUNA.ctd.level2 = fullfile(cfg.yeardir,'SBE16','CEO2_2017-2018_ctd033m_L2_v2.csv');% QC'd, averaged
    % Calibration file information - created BEFORE deployment
    cfg.SUNA.cal1.calfile = 'SNA0562J.cal';
    cfg.SUNA.cal1.caldate = '03-Aug-2017';
    cfg.SUNA.cal1.lampsec = 309390;    %cum_lampon_sec when (or very near when) cal file was created
    % Calibration file information - created AFTER recovery
    cfg.SUNA.cal2.calfile = 'SNA0562K.cal';
    cfg.SUNA.cal2.caldate = '10-Aug-2019';
    cfg.SUNA.cal2.lampsec = 466822;    %cum_lampon_sec when (or very near when) cal file was created
  case 'CEO_2018'
    % calibration information and files
    fprintf('dont know yet...\n')
    % ctd information
    cfg.SUNA.ctd.type  = 'SBE16';
    cfg.SUNA.ctd.SN    = cfg.mooring.SN.(cfg.SUNA.ctd.type);
    cfg.SUNA.ctd.level1  = fullfile(cfg.yeardir,'SBE16','CEO2_2018-2019_SBE16_4604.cnv');    % Full resoltuion - use to determine time lag
    cfg.SUNA.ctd.level2  = fullfile(cfg.yeardir,'SBE16','CEO2_2018-2019_ctd033m_L2_v2.csv'); % QC'd, averaged
    % calibration information and files
    cfg.SUNA.cal.postcal = struct();
    cfg.SUNA.cal.postcal.files  = 'D2019231_runwithDI.CSV';
    cfg.SUNA.cal.postcal.ctd_file   = '';  % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.postcal.NO3_uM     = 0.0;  % NO Nitrate standard available, ONLY DI measurements this year!
    % Calibration file information - created BEFORE deployment
    cfg.SUNA.cal1.calfile = 'SNA0801E.cal';
    cfg.SUNA.cal1.caldate = '24-Jul-2018';
    cfg.SUNA.cal1.lampsec = 416500;         %cum_lampon_sec when (or very near when) cal file was created
    if isfield(cfg.SUNA,'cal2')
      cfg.SUNA = rmfield(cfg.SUNA,'cal2');
    end
  case 'CEO_2019'
    % calibration information and files
    % predeployment measurements
    cfg.SUNA.cal.precal          = struct();
    cfg.SUNA.cal.precal.files    = {'D2019222.CSV'};
    cfg.SUNA.cal.precal.NO3_uM   = 29.967;  % NO Nitrate standard available, ONLY DI measurements this year!
    cfg.SUNA.cal.precal.ctd_file = '';      % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.precal.comment  = 'Acros Potassium nitrate, purity 99.999% (Lot A0345962). Final concentration NO3 standard = 29.967 然 in MQ'; % Email from Mette Kaufman <mrkaufman@alaska.edu> on Nov 3 2020
    % post deployment measurements
    cfg.SUNA.cal.postcal          = struct();
    cfg.SUNA.cal.postcal.files    = {'D2020285sn562_standards.CSV'};
    cfg.SUNA.cal.postcal.ctd_file = '';    % CTD data taken coincident with reference measurements
    cfg.SUNA.cal.postcal.NO3_uM   = 15.486;   % Nitrate standard concentration in uM. Set to NaN if unknown - defaults to 30uM 
    % Email correspondence with Ana & Meta | Subject: Nitrate standard
    % Question | Date: Nov 11 2020
    % In 2020, I gave Seth a 35.00  然  standard in MQ (that I made
    % onboard) for calibration during Spring, Summer cruises Two standards
    % were available during the fall cruise: The same standard as above,
    % and a 15.486 然  standard in MQ that Mette had made. Ethan borrowed
    % one, but did not tell Emily or Mette which one he grabbed. Ethan
    % might have been trying to calibrate the ISUS for the underway...not
    % sure For what I remember helping Seth with calibration, the SUNA can
    % only take MQ (to get a spectrum for 0 nitrate) and a high standard
    % (to get the spectrum for the high end range).


    % ctd information
    cfg.SUNA.ctd.type    = 'SBE16';
    cfg.SUNA.ctd.SN      = cfg.mooring.SN.(cfg.SUNA.ctd.type);
    cfg.SUNA.ctd.level1  = fullfile(cfg.yeardir,'SBE16','CEO2_2019-2020_SBE16_6670.cnv'); % Full resoltuion - use to determine time lag
    cfg.SUNA.ctd.level2  = fullfile(cfg.yeardir,'SBE16','CEO2_2019-2020_ctd033m_L2_v3.csv'); % QC'd, averaged
   % Calibration file information - created BEFORE deployment
    cfg.SUNA.cal1.calfile = 'SNA0562K.cal';
    cfg.SUNA.cal1.caldate = '10-Aug-2019';
    cfg.SUNA.cal1.lampsec = 466822;    %cum_lampon_sec when (or very near when) cal file was created
  otherwise
    fprintf('Project not set up yet: %s\n',cfg.project)
    keyboard
end

%% Read calibration data
caldir = fullfile(cfg.datadir,'CAL');
% Build file path
if exist(fullfile(caldir,cfg.SUNA.cal1.calfile),'file')
  calfile = fullfile(caldir,cfg.SUNA.cal1.calfile);
elseif exist(fullfile(caldir,'deployment',cfg.SUNA.cal1.calfile),'file')
  calfile = fullfile(caldir,'deployment',cfg.SUNA.cal1.calfile);
else
  fprintf('could not find %s...\n',cfg.SUNA.cal1.calfile)
  keyboard
end
[cfg.SUNA.cal1.caldata, ~] = SUNA_read_calfile(calfile);
% check if secondary calibration file exists
if isfield(cfg.SUNA,'cal2')
  % build file path
  if exist(fullfile(caldir,cfg.SUNA.cal2.calfile),'file')
    calfile = fullfile(caldir,cfg.SUNA.cal2.calfile);
  elseif exist(fullfile(caldir,'deployment',cfg.SUNA.cal2.calfile),'file')
    calfile = fullfile(caldir,'deployment',cfg.SUNA.cal2.calfile);
  else
    fprintf('could not find %s...\n',cfg.SUNA.cal2.calfile)
    keyboard
  end
  [cfg.SUNA.cal2.caldata, ~] = SUNA_read_calfile(calfile);
end

%% SUNA Nitrate V2 Instrument specifications
cfg.SUNA.model      = 'V2';
cfg.SUNA.pathlength = 10;
cfg.SUNA.spectral_range = [190 370]; % nm
cfg.SUNA.light_lifetime = 900;       % hours
cfg.SUNA.microM_to_mgL  = 0.014007;  % 1 然 nitrogen = 0.014007 mg/L nitrogen

%% Nitrate measurement accuracy (manual 1.4.1)
% From user manual
% User manual Submersible Ultraviolet Nitrate Analyzer SUNA V2, 2017
% Note: these are all specs for seawater!
cfg.SUNA.range.uM  = [0    1000; 1000 2000; 2000 3000]; % 然
cfg.SUNA.range.mgL = [0  14; 14 28; 28 42];  % mgN/L
switch cfg.SUNA.pathlength
  case 10 % mm
    cfg.SUNA.accuracy.best    = 2.0;          % [ 然 ] Best = 0.028 mgN/L (2 然)
    cfg.SUNA.accuracy.percent = [10; 15; 20]; % [ % ] i.e. for [range1; range2; range3]
  case 5  % mm
    cfg.SUNA.accuracy.best    = 4.0;          % [ 然 ] Best = 0.056 mgN/L (4 然)
    cfg.SUNA.accuracy.percent = [10; 15; 15]; % [ % ] i.e. for [range1; range2; range3]
  otherwise
    fprintf('pathlength not recognized, only 5 or 10mm possible\n')
    keyboard
end
%% Nitrate measurement precision (manual 1.4.2)
cfg.SUNA.prec = struct();
% Seawater or fresh water with T-S correction
cfg.SUNA.prec.TScor.shortterm      = 0.3; % [ 然 ] 0.3 然 (0.004 mgN/L)    = short-term precision (3 sigma) and limit of detection
cfg.SUNA.prec.TScor.drift_lamptime = 0.3; % [ 然 ] < 0.3 然 (< 0.004 mgN/L)= change ("drift") per hour of lamp time
cfg.SUNA.prec.TScor.limit_quant    = 1.0; % [ 然 ] 1.0 然 (0.014 mgN/L);   = limit of quantification
% Seawater (0-40 PSU)
cfg.SUNA.prec.SW.shortterm      = 2.4; % [ 然 ] 2.4 然 (0.034 mgN/L)       = short-term precision (3 sigma) and limit of detection
cfg.SUNA.prec.SW.drift_lamptime = 1.0; % [ 然 ] < 1.0 然 (< 0.014 mgN/L)   = change ("drift") per hour of lamp time
cfg.SUNA.prec.SW.limit_quant    = 8.0; % [ 然 ] 8.0 然 (0.112 mgN/L)       = limit of quantification

%% Data precision for writing data to files
cfg.SUNA.prec.fmt = '%.1f'; % 1 decimal place

%% NetCDF standard metadata
cfg.metadata.instrument_SUNA.core_variable      = 'dissolved_nutrients';
cfg.metadata.instrument_SUNA.parameter          = 'nitrate';
cfg.metadata.instrument_SUNA.unit               = 'umol/L';
cfg.metadata.instrument_SUNA.long_name          = 'Sea-Bird Scientific Submersible Ultraviolet Nitrate Analyser (SUNA) V2';
cfg.metadata.instrument_SUNA.instrument_type    = 'Spectrophotometer';
cfg.metadata.instrument_SUNA.pathlength_cm      = cfg.SUNA.pathlength;
cfg.metadata.instrument_SUNA.make_model         = ['SUNA ' cfg.SUNA.model];
cfg.metadata.instrument_SUNA.serial_number      = cfg.mooring.SN.SUNA;
cfg.metadata.instrument_SUNA.calibration_date   = cfg.SUNA.cal1.caldate;
cfg.metadata.instrument_SUNA.calibration_report = cfg.SUNA.cal1.calfile;
cfg.metadata.instrument_SUNA.accuracy           = cfg.SUNA.accuracy.best;    
cfg.metadata.instrument_SUNA.precision          = cfg.SUNA.prec.TScor.shortterm;  % Seawater or fresh water with T-S correction
cfg.metadata.instrument_SUNA.temperature_range  = [0 35];                         % temperature range, operation
cfg.metadata.instrument_SUNA.salinity_range     = [0 40];                         % salinity range

switch cfg.mooring.SN.SUNA  
  case '562' % CEO 2015-2016, 2017-2018, 2019-2020
    cfg.metadata.instrument_SUNA.factory_calibrated = '04-Apr-2015';
  case '801' % CEO 2016-2017, 2018-2019
    cfg.metadata.instrument_SUNA.factory_calibrated = '28-Jun-2016';
  case '840' % CEO 2020-2021
    cfg.metadata.instrument_SUNA.factory_calibrated = '06-Jun-2017';
end

%% metadata for instrument_ctd
cfg.metadata.instrument_ctd.long_name     = ['Sea-Bird Scientific ' cfg.SUNA.ctd.type];
cfg.metadata.instrument_ctd.make_model    = ['SBE ' cfg.SUNA.ctd.type];
cfg.metadata.instrument_ctd.serial_number = cfg.SUNA.ctd.SN;
if isfield(cfg.SUNA.ctd,'level2')
  cfg.metadata.instrument_ctd.filename    = cfg.SUNA.ctd.level2;
else
  cfg.metadata.instrument_ctd.filename    = cfg.SUNA.ctd.level1;
end

switch cfg.SUNA.ctd.type
  case 'SBE16'
    cfg.metadata.instrument_ctd.accuracy_conductivity = 0.005;
    cfg.metadata.instrument_ctd.accuracy_temperature  = 0.005;
  case 'SBE37'
end

end