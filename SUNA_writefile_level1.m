function SUNA_writefile_level1(cfg,data_raw)
%FUNCTION SUNA_writefile_level1 
%
%  Syntax:
%     SUNA_writefile_level1(cfg,data_raw)
%  
%  Description:
%     Writes level 1 data to CSV file containing ASCII version of full
%     resolution raw data from that deployment year. 
%     File naming convention: 
%     [mooring_code]_[year-range]_SUNAv2_Nitrate_L[processinglevel]_v[data version].csv.
%
%  References: 
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% 0 | Script set up
opt = struct();
opt.L1_version = 1; % Level1 data version - default to 1, change if file exists
opt.comment = '%';
opt.delm    = ',';

%% Generate filenames
level1_filename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L1_v' num2str(opt.L1_version) '.csv']);

while exist(level1_filename,'file')
  opt.L1_version = opt.L1_version + 1;
 level1_filename = fullfile(cfg.datadir,[cfg.metadata.global.mooring_code '_' num2str(cfg.year) '-' num2str(cfg.year+1) '_SUNAv2_Nitrate_L1_v' num2str(opt.L1_version) '.csv']);
end

%% Write Level1 file, ASCII version of raw data, full resolution
l1_choice = input(' Do you want to write Level1 data to file? <1/0> ');
if isempty(l1_choice); l1_choice = 0; end % default to NO
if l1_choice
  fprintf('Writing Level1 data to: %s\n',level1_filename);
  writetable(data_raw,level1_filename)
end
end %% function [cfg, t] = SeapHOx_find_zero(cfg,data)