function SUNA_write_calfile(cfg,cal_data)
% FUNCTION SUNA_write_calfile
%
%  Syntax:
%    SUNA_write_calfile(cfg,cal_data)
%  
%  Description:
%    Write SUNA .cal file from reference data read from CSV file
%
%  See also:
%    User manual Submersible Ultraviolet Nitrate Analyzer SUNA V2, 2017
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Check input arguments okay
if nargin < 2
  error('must pass through "cfg" structure and "cal_data" structure that contains data from refence measurements that will be used to write new calibration file')
end

%% 1 | Open file for writing
new_calfile = fullfile(cfg.datadir,'CAL',cal_data.calfile);
fprintf(' Writing new calibration file: %s\n',new_calfile)

fileID = fopen(new_calfile,'w');
if fileID < 0
  fprintf(' *** Error opening file %s\n',new_calfile)
  keyboard
end

%% 2 | Write header and reference data to file
% write header
fprintf(fileID,'%s\n',cal_data.header{:});  

% convert table to cell array to handle difference variable types
writedata = table2cell(cal_data.caldata); 
% transpose because fprintf function prints data columnwise
writedata = writedata'; 
% Define format for writing
fmt = '%s,%.2f,%.8f,%.8f,%.8f,%.8f\n';

fprintf(fileID,fmt,writedata{:});     % write data
fclose(fileID);   

end