function [cal, T_CAL] = SUNA_read_calfile(calfile)
% FUNCTION SUNA_read_calfile
%
%  Syntax:
%    cfg = SUNA_read_calfile(calfile)
%  
%  Description:
%    Read SUNA .cal file
%
%  See also:
%    User manual Submersible Ultraviolet Nitrate Analyzer SUNA V2, 2017
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 
if nargin < 1
  error('must pass through calibration file')
end
%% Get calibration file import options
calopts = detectImportOptions(calfile,'FileType','text');

%% Read calibration temperature from .cal header
fileid = fopen(calfile,'r');
for nl = 1:calopts.DataLine-1
  str = fgetl(fileid);
  str = strsplit(str,'H,');
  if contains(str{2},'T_CAL')
   T_CAL = str2double(erase(str{2},'T_CAL'));
  end
end
fclose(fileid);

%% Read T_CAL from .cal header
cal = readtable(calfile,calopts);

end