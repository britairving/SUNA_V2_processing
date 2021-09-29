function SUNA_processing_workflow
%FUNCTION SUNA_processing_workflow 
%
%  Syntax:
%    SUNA_processing_workflow
%  
%  Description:
%    This toolbox was developed to process data collected with a SUNA V2
%    deployed on a mooring.
%    
%    Processing scripts were written in MATLAB version R2017b using
%    Microsoft Windows 10 Pro but should be portable to other MATLAB
%    versions and operating systems. The toolbox will be published to
%    Github in the future.
%  
%    data_raw | table containing raw data
%    meta_raw | structure containing metadata for each variable in data_raw
%    data_pre | structure containing selected SUNA and CTD data variables with standardized names and units. 
%    meta_pre | structure containing source, and conversion if performed on each variable in data_pre.
%    data_proc | structure containing full resolution processed data
%    meta_proc | structure containing units and information on processing steps performed on each variable in data_proc.
%    data_avg | structure containing burst resolution processed data 
%    meta_avg | structure containing units and information on processing steps performed on each variable in data_avg.
%
%  References: 
%    https://github.alaska.edu/bkirving/SUNA_Nitrate_processing
%
%  Dependencies:
%    https://github.com/britairving/data_analysis_tools
%    TEOS-10 GSW Toolbox: https://github.com/TEOS-10/GSW-Matlab
%    M_MAP Toolbox: https://www.eoas.ubc.ca/%7Erich/map.html
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Define basic paths
cfg.year    = 2017;
cfg.project = ['CEO_' num2str(cfg.year)];
if ispc
  cfg.yeardir = fullfile('F','CEO','DataByDeploymentYear', num2str(cfg.year));
elseif ismac
  cfg.yeardir = fullfile(filesep,'Volumes','CEO','CEO','DataByDeploymentYear', num2str(cfg.year));
end
cfg.datadir = fullfile(cfg.yeardir,'SUNA',filesep);
cfg.save_figures = true; % Saves all figures to cfg.datadir;

%% 1 | Mooring configuration
cfg = mooring_config(cfg);

%% 2 | Processing configuration
% Build initial metadata by outlining project details, what variables,
% calibration files, processing steps, etc etc
cfg = SUNA_config(cfg);

%% 3 | Process data following best practices, or load already processed data
if exist(cfg.SUNA.path.data_proc,'file')
  fprintf(' Loading processed data to %s\n',cfg.SUNA.path.data_proc)
  % Only load the data, not the cfg structure
  load(cfg.SUNA.path.data_proc,'data_raw','meta_raw','data_pre','meta_pre','data_proc','meta_proc','data_avg','meta_avg')
else
  %% 4 | Read in raw data
  % Read in all datafiles from directory chosen in processing options
  if exist(cfg.SUNA.path.data_raw,'file')
    fprintf('Loading file: %s\n',cfg.SUNA.path.data_raw);
    load(cfg.SUNA.path.data_raw);
    % update cfg structure
    cfg = CEO_config(cfg);
    cfg = SUNA_config(cfg);
  else
    [cfg,data_raw,meta_raw] = SUNA_read_data(cfg);
    % Write L1 datafile
    SUNA_writefile_level1(cfg,data_raw);
  end
  
  %% 5 | Preprocess data
  % Convert raw instrument variable names to user friendly standardized names
  [cfg,data_pre,meta_pre] = SUNA_preprocessing(cfg,data_raw,meta_raw);
  
  %% 6 | Process data
  [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_processing(cfg,data_pre,meta_pre);
  fprintf(' Saving processed data to %s\n',cfg.SUNA.path.data_proc)
  save(cfg.SUNA.path.data_proc,'data_raw','meta_raw','data_pre','meta_pre','data_proc','meta_proc','data_avg','meta_avg','cfg','-v7.3');
end

%% 7 | Correct data
% Process data following best practices
if exist(cfg.SUNA.path.data_cor,'file')
  fprintf(' Loading corrected data to %s\n',cfg.SUNA.path.data_cor)
  load(cfg.SUNA.path.data_cor);
else
  [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_calibration_offset_drift_correction(cfg,data_proc,meta_proc,data_avg,meta_avg);
  fprintf(' Saving corrected data to %s\n',cfg.SUNA.path.data_cor)
  save(cfg.SUNA.path.data_cor,'data_raw','meta_raw','data_pre','meta_pre','data_proc','meta_proc','data_avg','meta_avg','cfg','-v7.3');
  % Write Level 2 and Level 3 files
  SUNA_writefile_level2_level3(cfg,data_proc,meta_proc,data_avg)
end

%% 8 | Uncertainty  
% TO DO....
% [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_error_propagation(cfg,data_proc,meta_proc,data_avg,meta_avg);
% fprintf(' Saving corrected data to %s\n',cfg.SUNA.path.data_cor)
% save(cfg.SUNA.path.data_cor,'data_raw','meta_raw','data_pre','meta_pre','data_proc','meta_proc','data_avg','meta_avg','cfg','-v7.3');

%% 9 | Plot Spectragram
SUNA_plot_data(data_proc,meta_proc,cfg);

%% 10 | QAQC
% Automated quality assurance and quality control tests following best
% practices (where available)
try
  [data_avg, cfg] = manual_timeSeries_QC(cfg,data_avg,{'NO3_uM_cor' 'NO3_uM' 'cdom' 'density'});
  fprintf(' Saving corrected data to %s\n',cfg.SUNA.path.data_qc)
  save(cfg.SUNA.path.data_qc,'data_raw','meta_raw','data_pre','meta_pre','data_proc','meta_proc','data_avg','meta_avg','cfg','-v7.3');
  % Write Level 4 - TO DO ( copy level3 file, essentially)
  %  - Data QC'd from expert
	%	 - Drift correction updated with discrete bottle measurements
  % SUNA_writefile_level4(cfg,data_proc,meta_proc,data_avg)
catch
  fprintf('Cannot find script... download manual_timeSeries_QC.m from https://github.com/britairving/data_analysis_tools\n')
end
end