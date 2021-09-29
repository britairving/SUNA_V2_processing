function SUNA_writefile_level2_level3(cfg,data_proc,meta_proc,data_avg)
%FUNCTION SUNA_writefile_level2_level3
%
%  Syntax:
%     SUNA_writefile_level2_level3(cfg,data_raw, meta_raw,data_proc,meta_proc,data_avg)
%  
%  Description:       
%    Writes level 2 data to CSV file containing ASCII version of full
%    resolution processed and corrected data.
%    Writes level 3 data to CSV file containing ASCII and NetCDF version of
%    burst resolution processed and corrected data.
% 
%    Level 2: ASCII version of processed data, full resolution
% 				-	upload into Researh Workspace  with metadata
% 				- Automatic QC tests
% 				- Nitrate corrected using temperature compensated salinity subtracted (TCSS) method of Sakamoto et al., 2009
% 				- Calibration drift corrected with available reference measurements and discrete bottle data
%    Level 3: ASCII and NetCDF version of smoothed data, burst resolution
%         - Updated to Level 2
% 				- Data averaged by taking mean of each segment, or burst. Usually 1 burst(multiple measurements)/hour. 
% 			  -	Flagged data are ignored in averaging.
% 				- Data smoothed with moving median over each window of 35 hours (e.g. approximately 35 burts)
%
%     File naming convention: 
%     [mooring_code]_[year-range]_SUNAv2_Nitrate_L[processinglevel]_v[data version].csv.
%
%  References: 
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Script set up
opt = struct();
opt.comment = '%';
opt.delm    = ',';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
opt.L2_version   = 1; % Level2 ASCII version - default to 1, change if file exists
opt.L3_version   = 1; % Level3 ASCII version - default to 1, change if file exists
opt.L3_ncversion = 1; % Level3 NetCDF data version - default to 1, change if file exists
opt.L3_fillval = -999;
opt.time_fmt    = 'yyyy-mm-ddTHH:MM:SS'; % ISO 8601:2004 extended date format
opt.time_fmt_nc = 'yyyymmddHHMMSS';
cfg.metadata.global.date_created = datestr(now,opt.time_fmt);
%% Generate filenames
level2_filename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L2_v' num2str(opt.L2_version) '.csv']);
while exist(level2_filename,'file')
  opt.L2_version = opt.L2_version + 1;
 level2_filename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L2_v' num2str(opt.L2_version) '.csv']);
end
% ASCII version
level3_filename   = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L3_v' num2str(opt.L3_version) '.csv']);
level3_ncfilename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L3_v' num2str(opt.L3_version) '.nc']);
while exist(level3_filename,'file')
  opt.L3_version = opt.L3_version + 1;
 level3_filename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L3_v' num2str(opt.L3_version) '.csv']);
end
% create netcdf filename
while exist(level3_ncfilename,'file')
  opt.L3_ncversion = opt.L3_ncversion + 1;
 level3_ncfilename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L3_v' num2str(opt.L3_ncversion) '.nc']);
end
%% Write Level2 file, ASCII version of processed data, full resolution

data_proc.DateTime = cellstr(datestr(data_proc.datenum,opt.time_fmt));
data_avg.DateTime  = cellstr(datestr(data_avg.datenum,opt.time_fmt));
% ALL variables that want Level2 file to have
vars      = struct();
vars.name = {'DateTime' 'NO3_uM'      'NO3_uM_cor'        'flag'         'dark'         'cum_lampon_sec'     'burst_index' 'pressure' 'temperature' 'salinity'};
vars.desc = {'DateTime' 'Nitrate_raw' 'Nitrate_corrected' 'Nitrate_flag' 'dark_counts'  'Cumulative_lamp_on' 'burst_index' 'Pressure' 'Temperature' 'Salinity'};
vars.ctd  = {'pressure' 'temperature' 'salinity' 'cdom' 'fluor_mgm3' 'turbidity' 'PAR' 'density' 'depth'};
dat      = struct();
datavg   = struct();
hdr      = {};
flds     = {};
units    = {};
for nv = 1:numel(vars.name)
  vname = vars.name{nv}; % variable name in data
  sname = vars.desc{nv}; % more descriptive name
  if isfield(data_proc,vname)
    % Fill data structure and variable header names
    hdr   = [hdr sname];
    flds  = [flds vname]; % keep track of which fields are available
    dat.(vname)    = data_proc.(vname);
    datavg.(vname) = data_avg.(vname);
    % Update units
    if strcmp(vname,'DateTime')
      units = [units 'yyyy-mm-ddTHH:MM:SS'];
    elseif strcmp(vname,'flag')
      units = [units 'nodim'];
    elseif contains(vname,'NO3')
      units = [units 'umol/L']; % micromolar
    else
      try
        units = [units strrep(meta_proc.(vname).units,' ','_')];
      catch
        fprintf('units not working\n')
        keyboard
      end
    end
  end
end

%% GENERATE FORMAT
fmt = cell(size(flds'));
for nvar = 1:numel(flds)
  field = flds{nvar};
  % Determine if at end of line
  if nvar == numel(flds)
    delm = '\n';
  else
    delm = opt.delm;
  end
  % change any character arrays to cell arrays
  if ischar(dat.(field))
    dat.(field) = cellstr(dat.(field));
  end
  % generate fmt
  if iscell(dat.(field)) || isdatetime(dat.(field))
    fmt{nvar} = ['%s' delm];
  elseif strcmp(field,'burst_index')
    fmt{nvar} = ['%.1f' delm];
  elseif contains(field,'NO3')
    fmt{nvar} = [cfg.SUNA.prec.fmt delm];
  elseif floor(dat.(field)(2)) == dat.(field)(2)
    fmt{nvar} = ['%d' delm];
  else
    fmt{nvar} = ['%.4f' delm];
  end
end
fmthdr = ['%s' opt.delm];
% 
%% Write Level2 and Level3 CSV files
%Level2 = ASCII version of processed data, full resolution
%Level3 = ASCII version of processed data, burst resolution
for nlev = 2:3
  if nlev == 2
    fprintf('Writing data to %s\n',level2_filename);
    fid = fopen(level2_filename,'w');
  elseif nlev == 3
    fprintf('Writing data to %s\n',level3_filename);
    fid = fopen(level3_filename,'w');
  end
  fprintf(fid,'%s\n',[opt.comment 'File Created: ' datestr(now)]);
  fprintf(fid,'%s\n',[opt.comment 'Mooring name: ' cfg.metadata.global.mooring_name]);
  fprintf(fid,'%s\n',[opt.comment 'Mooring code: ' cfg.metadata.global.mooring_code]);
  fprintf(fid,'%s\n',[opt.comment 'Latitude: '     num2str(cfg.mooring.latitude)   ' [decimal degrees N]']);
  fprintf(fid,'%s\n',[opt.comment 'Longitude: '    num2str(cfg.mooring.longitude) ' [decimal degrees W]']);
  fprintf(fid,'%s\n',[opt.comment 'Deployed: '     datestr(cfg.mooring.deploydate,'yyyy-mm-ddTHH:MM:SS')]);
  fprintf(fid,'%s\n',[opt.comment 'Recovered: '    datestr(cfg.mooring.recoverdate,'yyyy-mm-ddTHH:MM:SS')]);
  fprintf(fid,'%s\n',[opt.comment 'SUNA Calibration file: ' meta_proc.NO3_uM_cor.calfile]);
  fprintf(fid,'%s\n',[opt.comment 'CTD file: ' cfg.metadata.instrument_ctd.filename]);
  for nf = 1:numel(hdr)
    if nf == numel(hdr)
      fprintf(fid,'%s\n',hdr{nf});
    else
      fprintf(fid,fmthdr,hdr{nf});
    end
  end
  for nf = 1:numel(units)
    if nf == numel(units)
      fprintf(fid,'%s\n',units{nf});
    else
      fprintf(fid,fmthdr,units{nf});
    end
  end
  %% LEVEL 2
  if nlev == 2 
    for nline = 1:numel(dat.DateTime)
      for nvar = 1:numel(flds)
        field = flds{nvar};
        if iscell(dat.(field))
          fprintf(fid,fmt{nvar},dat.(field){nline});
        else
          fprintf(fid,fmt{nvar},dat.(field)(nline));
        end
      end
    end
  %% LEVEL 3
  elseif nlev == 3 
    for nline = 1:numel(datavg.DateTime)
      for nvar = 1:numel(flds)
        field = flds{nvar};
        if iscell(datavg.(field))
          fprintf(fid,fmt{nvar},datavg.(field){nline});
        else
          fprintf(fid,fmt{nvar},datavg.(field)(nline));
        end
      end
    end
  end
  % Close file
  fclose(fid);
end % Write Level2 file, then write level3 file

%% Write Level3 NetCDF file
fprintf('Writing NetCDF file: %s\n',level3_ncfilename);
% reformat time as integer
datavg.DateTime = cellstr(datestr(data_avg.datenum,opt.time_fmt_nc));
datavg.DateTime = str2double(datavg.DateTime);
try
  % create new file for writing
  scope = netcdf.create(level3_ncfilename,'NETCDF4');
  % Define dimensions:
  dimidT = netcdf.defDim(scope,'time',numel(datavg.DateTime));
  % Write global attributes
  globid = netcdf.getConstant('GLOBAL');
  global_fields = fieldnames(cfg.metadata.global);
  for ng = 1:numel(global_fields)
    netcdf.putAtt(scope,globid,global_fields{ng},cfg.metadata.global.(global_fields{ng}));
  end
  % Write SUNA instrument specs
  varid = netcdf.defVar(scope,'instrument_SUNA','char',[]);
  inst_fields = fieldnames(cfg.metadata.instrument_SUNA);
  for ni = 1:numel(inst_fields)
    netcdf.putAtt(scope,varid,inst_fields{ni},cfg.metadata.instrument_SUNA.(inst_fields{ni}));
  end
  % writeinstrument_ctd specs
  varid = netcdf.defVar(scope,'instrument_ctd','char',[]);
  ctd_fields = fieldnames(cfg.metadata.instrument_ctd);
  for ni = 1:numel(ctd_fields)
    netcdf.putAtt(scope,varid,ctd_fields{ni},cfg.metadata.instrument_ctd.(ctd_fields{ni}));
  end
catch
  fprintf('Could not write global metadata or instrument info to netcdf file\n')
  keyboard
end
try
  % Loop through variables
  for nvar = 1:numel(flds)
    field = flds{nvar};
    if strcmp(field,'DateTime')
      
      varid = netcdf.defVar(scope,hdr{nvar},'double',dimidT);
      netcdf.putVar(scope,varid,datavg.(field));
      netcdf.putAtt(scope,varid,'instrument','instrument_SUNA');
      netcdf.putAtt(scope,varid,'long_name' ,'Date time integer');
      netcdf.putAtt(scope,varid,'units',opt.time_fmt_nc);
      % also write MATLAB datenum
      varid = netcdf.defVar(scope,'matlab_datenum','double',dimidT);
      netcdf.putVar(scope,varid,data_avg.datenum);
      netcdf.putAtt(scope,varid,'instrument','instrument_SUNA');
      netcdf.putAtt(scope,varid,'long_name' ,'MATLAB serial date number');
      netcdf.putAtt(scope,varid,'units'     ,meta_proc.datenum.units);
    else
      % get variable ID
      varid = netcdf.defVar(scope,hdr{nvar},'double',dimidT);
      % Define fill value and replace NaNs with fill value
      datavg.(field)(isnan(datavg.(field))) = opt.L3_fillval;
      netcdf.defVarFill(scope,varid,false,opt.L3_fillval);
      
      %if strcmp(field,'temperature')
      %  keyboard
      %end
      % write data to netcdf file
      netcdf.putVar(scope,varid,datavg.(field));

      % write data attributes to netcdf file
      if ismember(field,vars.ctd)
        netcdf.putAtt(scope,varid,'instrument','instrument_ctd');
      else
        netcdf.putAtt(scope,varid,'instrument','instrument_SUNA');
      end
      if isfield(meta_proc.(field),'standard_name')
        netcdf.putAtt(scope,varid,'standard_name',meta_proc.(field).standard_name);
      end
      if isfield(meta_proc.(field),'meaning') % meaning is usually more descriptive
        netcdf.putAtt(scope,varid,'long_name',meta_proc.(field).meaning);
      elseif isfield(meta_proc.(field),'long_name')
        netcdf.putAtt(scope,varid,'long_name',meta_proc.(field).long_name);
      end
      if isfield(meta_proc.(field),'units')
        if strcmp(meta_proc.(field).units,'uM')
          meta_proc.(field).units = 'umol/L';
        end
        netcdf.putAtt(scope,varid,'units',meta_proc.(field).units);
      end
      if strcmp(cfg.SUNA.ctd.type,'SBE37') && (contains(field,'temperature','IgnoreCase',1) || contains(field,'salinity','IgnoreCase',1)) 
        fprintf('write cal info\n')
        keyboard
      end
      % write flag meanings
      if contains(field,'flag')
        netcdf.putAtt(scope,varid,'flag_meanings',cfg.metadata.flag_meanings);
        netcdf.putAtt(scope,varid,'flag_values',  cfg.metadata.flag_values);
      end
      % write as much metadata for corrected nitrate as possible 
      if contains(field,'cor')
        %fprintf('updated corrected fields  -- references, etc\n')
        refs = '';
        for ncal = 1:numel(meta_proc.(field).offsets.type)
          if contains(meta_proc.(field).offsets.type{ncal},'discrete')
            ref_string = ['Calibration point ' num2str(ncal) ' = Discrete sample(s) from ' meta_proc.(field).offsets.citation{contains(meta_proc.(field).offsets.type,'discrete')}];
          elseif contains(meta_proc.(field).offsets.type{ncal},'post-recovery')
            ref_string = ['Calibration point ' num2str(ncal) ' = absolute offset from updated reference postrecovery.cal'];
          elseif contains(meta_proc.(field).offsets.type{ncal},'pre-deployment')
            ref_string = ['Calibration point ' num2str(ncal) ' = absolute offset from interpolated pre-deployment reference from ' meta_proc.(field).reference{ncal}];
          elseif  contains(meta_proc.(field).offsets.type{ncal},'zero')
            ref_string = ['Calibration point ' num2str(ncal) ' = zero correction ' meta_proc.(field).offsets.reference{ncal}];
          else
            fprintf('not sure what to do in this case..\n')
            keyboard
          end
          if ncal == 1
            refs = ref_string;
          else
            refs = [refs ', ' ref_string];
          end
          if ~contains(refs,'Sakamoto')
            refs = ['Sakamoto et al. (2009) Temperature compensated salinity subtracted (TCSS) doi.org/10.4319/lom.2009.7.132 & Randelhoff et al. (2016) doi.org/10.1002/2016JC011779, ' refs];
          end
        end
        netcdf.putAtt(scope,varid,'references',refs);
      end % Corrected data
    end % IF Datetime
  end % Loop through variables to write to NetCDF file
catch
  fprintf('could not write file...%s\n',field)
  keyboard
end

netcdf.close(scope)
ncdisp(level3_ncfilename)
keyboard
end %% function [cfg, t] = SeapHOx_find_zero(cfg,data)