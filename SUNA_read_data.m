function [cfg,data,meta] = SUNA_read_data(cfg,file)
% FUNCTION SUNA_read_data
%
%  Syntax:
%    cfg = SUNA_read_data(cfg)
%    "cfg" is a structure containing information about the deployment and
%    the instrument.
%  
%  Description:
%    Reads SUNA csv files in the SUNA_L0 directory. 
%
%  See also:
%    User manual Submersible Ultraviolet Nitrate Analyzer SUNA V2, 2017
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% file was passed through
if nargin > 1
  read_single_file = 1;
else
  read_single_file = 0;
end

%% Define CSV column fieldnames
% from the user manual section 6.1.2.1 "File Types"
cols_spec  = strcat('spec',compose('%d',1:256)); % generate fieldnames of spectrum channels
cols_start = {'hdr' 'yyyydoy' 'hr' 'NO3_uM' 'NO3_mgNL' 'abs_254nm' 'abs_350nm' 'bromide_mgL' 'spec_avg' 'dark' 'int_fact'};
cols_end   = {'int_temp' 'spec_temp' 'lamp_temp' 'cum_lampon_sec' 'rel_hum' 'volt' 'lamp_volt' 'int_volt' 'current' 'fit_aux1' 'fit_aux2' 'fit_base1' 'fit_base2' 'fit_RMSE'...
              'ctd_time' 'ctd_sal' 'ctd_temp' 'ctd_press' 'check_sum'};
cols = [cols_start cols_spec cols_end];
% add units/description
cols_start_units = {'Header and serial number' 'yearDOY' 'hours' 'uM' 'mgN/L' 'nm' 'nm' 'mg/L' 'spectrum average' 'dark value used for fit' 'integration time factor'};
cols_spec_units  = strcat('spectrum channel',compose('%d',1:256)); % generate fieldnames of spectrum channels
cols_end_units   = {'degC' 'degC' 'degC' 'sec' '%' 'V' 'V' 'V' 'mA' 'fit aux 1' 'fit aux 2' 'fit base 1' 'fit base 2' 'fit RMSE' 'sec' 'PSU' 'degC' 'dbar' 'nodim'};
cols_units       = [cols_start_units cols_spec_units cols_end_units];

%% Check for data directory existance
if read_single_file == 0
  raw_dir = fullfile(cfg.datadir,'SUNA_L0');
  if ~exist(raw_dir,'dir')
    fprintf('Directory not found: %s\n',raw_dir)
    keyboard
  end
  
  %% Get mooring deployment date in yyyyDDD format (year day-of-year)
  if isfield(cfg,'mooring')
    deploydate = datetime( cfg.mooring.deploydate,'Format','yyyyDDD');
  else % default to today
    deploydate = datetime('now','Format','yyyyDDD');
  end
  if isfield(cfg,'mooring') && isfield(cfg.mooring,'recoverdate')
    recoverdate = datetime( cfg.mooring.recoverdate,'Format','yyyyDDD');
  else % default to today
    recoverdate = datetime('now','Format','yyyyDDD');
  end
  % convert to numeric value instead of datetime
  deploydate  = str2double(char(deploydate));
  recoverdate = str2double(char(recoverdate));
  
  %% Read all raw files
  files = dir(fullfile(raw_dir,'*.CSV'));
  % remove files with erroneous names
  ignore_files = contains({files.name},{'..' '._'});
  files(ignore_files)       = [];
  files([files.bytes] == 0) = []; % throw out empty files
  % skip files containing data before deployment
  yyyyDDD = str2double(erase({files.name},{'D' '.CSV'}));
  files(yyyyDDD < deploydate | yyyyDDD > recoverdate)  = [];
else
  files = struct();
  for nf = 1:size(file,1)
    [folder,name,ext] = fileparts(file(nf,:));
    files(nf).folder = folder;
    files(nf).name   = [name ext];
  end
end

%% Initialize meta structure 
% stores filenames, pertinent information
meta = struct();
% add column names and units to meta 
for nc = 1:numel(cols)
  meta.(cols{nc}).units = cols_units{nc};
end

% save files to meta structure
meta.csv_filenames = {files.name}';
for nf = 1:numel(files)
  filenam = fullfile(files(nf).folder,files(nf).name);
  fprintf('Reading raw SUNA datafile: %s\n',files(nf).name)
  opts = detectImportOptions(filenam);
  if numel(opts.VariableNames) == 286
    opts.VariableNames = cols;
  else
    fprintf('... unexpected first column (datetime) exists in file...\n')
    fprintf('... may cause erros ...\n')
    opts.VariableNames = ['datetime' cols];
  end
  % read in CTD variables as NaNs instead of ''
  opts.VariableTypes(contains(opts.VariableNames,{'ctd_sal' 'ctd_temp' 'ctd_press'})) = {'double'};
  %% First read header informatin
  inst  = struct(); % initialize inst and store all coefficients
  fileid = fopen(filenam,'r');
  for nl = 1:opts.DataLine-1
    try
      str = fgetl(fileid);
      str = strsplit(str,',');
      hdr_field = strrep(str{2},'-','_');
      hdr_field = regexprep(hdr_field,' +',' '); % replace multiple spaces with single space
      hdr_field = strtrim(hdr_field);
      hdr_field = strrep(hdr_field,' ','_');
      inst.(hdr_field) = strtrim(str{3});
    catch
      fprintf('error here\n')
      keyboard
    end
  end
  fclose(fileid);

  %% Next read data
	d = readtable(filenam,opts);
  d.filename = repmat(files(nf).name,size(d.hdr));
  if ~read_single_file
    try
      if d.yyyydoy(end) < deploydate || d.yyyydoy(end) > recoverdate
        % remove fieldname
        meta.csv_filenames(contains(meta.csv_filenames,files(nf).name)) = [];
        continue
      end
    catch
      fprintf('stopped here: %s\n',filenam);
      keyboard
    end
    % check to make sure header information is the same
    if exist('data','var') % i.e. one file has already been read
      hdr_fields = fieldnames(inst);
      for ncal = 1:numel(hdr_fields)
        hdr_field = hdr_fields{ncal};
        if ~strcmp(inst_prev.(hdr_field),inst.(hdr_field))
          fprintf('Header info different than previous! %s\n',filenam)
          fprintf('  previous file %s = %s\n',hdr_field,inst_prev.(hdr_field))
          fprintf('  this file     %s = %s\n',hdr_field,inst.(hdr_field))
          skip_file = input('Do you want to skip this file? <1/0> ');
          if isempty(skip_file); skip_file = 1; end
          if skip_file
            meta.csv_filenames(contains(meta.csv_filenames,files(nf).name)) = [];
            inst = inst_prev;
            continue
          else
            fprintf('stopped here...\n')
            keyboard
          end
        end
      end
    end
  end
  inst_prev = inst; % update so can check next file has the same header information
  if ~exist('data','var')
    data = d;
  else
    data  = [data;d];
  end
end
% Add datenum and datetime
data.datetime = datetime(string(data.yyyydoy),'InputFormat','yyyyDDD') + hours(data.hr);
data.datenum  = datenum(data.datetime);
% store CSV header information to meta structure
meta.csv_headers = inst;

% if ~strcmp(meta.csv_headers.Calibration_File,cfg.SUNA.cal1.calfile)
%   fprintf('Different calibration files...\n')
%   keyboard
% end

%% Save data
if ~read_single_file
  meta_raw = meta;
  data_raw = data;
  fprintf('Calibration file for %s = %s\n',cfg.project,meta_raw.csv_headers.Calibration_File)
  fprintf('Saving raw data to file: %s\n',cfg.SUNA.path.data_raw)
  save(cfg.SUNA.path.data_raw,'data_raw','meta_raw','cfg','-v7.3');
end

end