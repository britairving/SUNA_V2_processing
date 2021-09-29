function ctd = read_ctd_profiles_by_project(cfg)
% FUNCTION read_ctd_profiles_by_project
%
%  Syntax:
%    ctd = read_ctd_profiles_by_project(project_name)
%  
%  Description:
%    Load/read in relevant CTD profiles for each cruise. Cruise, Cast, and
%    Station must all be the same as those equivalent in cfg.path.calfile.
%
%  Inputs:
%    cfg | structure containing at least cfg.project, and cfg.paths to
%          point to where relevant ctd files are. 
%
%  Outputs:
%    ctd | structure containing CTD profile data.
%       ctd.Cruise   = cell array with cruise name
%       ctd.Cast     = double array of cast numbers
%       ctd.dnum     = double array with MATLAB's datenum
%       ctd.Station  = cell array with Station name
%       ctd.density  = double array of density  [kg/m^3]
%       ctd.pressure = double array of pressure [ dbar ]
%       ctd.salinity = double array of salinity [ psu  ]
% 
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
% fprintf('STANDARDIZE REQUIRED CTD FORMAT -- SO LOADING IS EASIER!\n')
% if ~strcmp(datestr(now,'yyyymmdd'),'20201117')
%   keyboard
% end
%% Read CTD profiles for specific project 
% Formatting greatly varies..... so user must figure this out
switch cfg.project
  case {'CEO_2016' 'CEO_2016_2017'}
    % AMBON2017
    ctd_file = fullfile(cfg.path.calcasts,'AMBON2015_2017','AMBON2017_ctd_L3_v1.csv');
    ctd1 = readtable(ctd_file);
    ctd_match = ismember(ctd1.Station,unique(calcasts.StationID));  %ctd_match = strcmp(ctd.Station,cfg.discrete_ref.StationID{1});
    ctd1 = ctd1(ctd_match,:);
    % calculate density, datenum, and datetime
    ctd1.density = sw_dens(ctd1.salinity__psu_,ctd1.temperature__C_,ctd1.pressure__dbar_);
    ctd1.dtime = datetime(strrep(ctd1.Date_Time,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:SS');
    ctd1.dnum  = datenum(ctd1.dtime);
    ctd1.pressure = ctd1.pressure__dbar_;
    ctd1.salinity = ctd1.salinity__psu_;
    % SKQ201612S
    ctd2 = cnv2mat_bi(fullfile(cfg.path.calcasts,'2016_CastC045')); % pass through directory, not filename
    ctd2.Cruise  = cellstr(repmat('SKQ201612S',size(ctd2.datenum)));
    ctd2.Cast    = repmat(45,size(ctd2.datenum));
    ctd2.Station = repmat(ctd2.station,size(ctd2.datenum));
    ctd2.density = sw_dens(ctd2.salinity,ctd2.temperature,ctd2.pressure);
    % add together
    ctd = struct();
    ctd.Cruise   = [ctd1.Cruise;   ctd2.Cruise];
    ctd.Cast     = [ctd1.Cast;     ctd2.Cast];
    ctd.dnum     = [ctd1.dnum;     ctd2.datenum];
    ctd.Station  = [ctd1.Station;  ctd2.Station];
    ctd.density  = [ctd1.density;  ctd2.density];
    ctd.pressure = [ctd1.pressure; ctd2.pressure];
    ctd.salinity = [ctd1.salinity; ctd2.salinity];
  case 'CEO_2015'
    % AMBON2015
    ctd_file = fullfile(cfg.path.calcasts,'AMBON2015_2017','AMBON2015_ctd_L3_v1.csv');
    ctd = readtable(ctd_file);
    ctd_match = ismember(ctd.Station,unique(calcasts.StationID));  %ctd_match = strcmp(ctd.Station,cfg.discrete_ref.StationID{1});
    ctd = ctd(ctd_match,:);
    % calculate density, datenum, and datetime
    ctd.density = sw_dens(ctd.salinity__psu_,ctd.temperature__C_,ctd.pressure__dbar_);
    ctd.dtime = datetime(strrep(ctd.Date_Time,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:SS');
    ctd.dnum  = datenum(ctd.dtime);
    ctd.pressure = ctd.pressure__dbar_;
    ctd.salinity = ctd.salinity__psu_;
  case {'CEO_2017' 'CEO_2017_2018'}
    % HLY1702 [https://web.whoi.edu/healy-2017/]
    ctd_file1 = fullfile(cfg.path.calcasts,'HLY1702','Calibrated_ctd_files','HLY1702_126.dcc');
    ctd1 = readtable(ctd_file1,'FileType','text');
    ctd1.Cast = repmat(126,size(ctd1.Pres));
    ctd1.Station = cellstr(repmat('S1',size(ctd1.Pres)));
    ctd1.dnum    = repmat(datenum(2017,09,09,23,36,00),size(ctd1.Pres));
    
    ctd_file2 = fullfile(cfg.path.calcasts,'HLY1702','Calibrated_ctd_files','HLY1702_127.dcc');
    ctd2 = readtable(ctd_file2,'FileType','text');
    ctd2.Cast = repmat(127,size(ctd2.Pres));
    ctd2.Station = cellstr(repmat('CEO',size(ctd2.Pres)));
    ctd2.dnum    = repmat(datenum(2017,09,10,01,12,00),size(ctd2.Pres));
    % merge both casts into single table
    ctd_HLY1702 = [ctd1; ctd2];
    % calculate density
    ctd_HLY1702.density  = sw_dens(ctd_HLY1702.Sal_1_,ctd_HLY1702.T90_1_,ctd_HLY1702.Pres);
    ctd_HLY1702.pressure = ctd_HLY1702.Pres;
    ctd_HLY1702.salinity = ctd_HLY1702.Sal_1_;
    ctd_HLY1702.Cruise   = cellstr(repmat('HLY1702',size(ctd_HLY1702.pressure)));
    % Read AMBON 2017 data
    ctd_file = fullfile(cfg.path.calcasts,'AMBON2015_2017','AMBON2017_ctd_L3_v1.csv');
    ctd3 = readtable(ctd_file);
    ctd_match = ismember(ctd3.Station,unique(calcasts.StationID));  %ctd_match = strcmp(ctd3.Station,cfg.discrete_ref.StationID{1});
    ctd3 = ctd3(ctd_match,:);
    % calculate density, datenum, and datetime
    ctd3.density = sw_dens(ctd3.salinity__psu_,ctd3.temperature__C_,ctd3.pressure__dbar_);
    ctd3.dtime = datetime(strrep(ctd3.Date_Time,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:SS');
    ctd3.dnum  = datenum(ctd3.dtime);
    ctd3.pressure = ctd3.pressure__dbar_;
    ctd3.salinity = ctd3.salinity__psu_;
    % create one structure
    ctd = struct();
    ctd.Cruise   = [ctd_HLY1702.Cruise;  ctd3.Cruise];
    ctd.Cast     = [ctd_HLY1702.Cast;    ctd3.Cast];
    ctd.dnum     = [ctd_HLY1702.dnum;    ctd3.dnum];
    ctd.Station  = [ctd_HLY1702.Station; ctd3.Station];
    ctd.density  = [ctd_HLY1702.density; ctd3.density];
    ctd.pressure = [ctd_HLY1702.pressure;ctd3.pressure];
    ctd.salinity = [ctd_HLY1702.salinity;ctd3.salinity];
    % SKQ201818S -- no ctd file but can add basic info from excel file
    % D:\HauriLab\CEOmooring\CalibrationSamples\original files\SKQ2018_18S_CEO_2018_deployment_calibration.xlsx
    t = [-1.5974	-1.572 -1.1799	2.122	1.9901];
    s = [32.6865	32.6726 32.3311	29.473	29.3767];
    d = [41.9 33.3	20.2	10.1	3.4];% depth
    p = gsw_p_from_z(-d,71.603);
    d = sw_dens(s,t,p);
    ctd.Cruise   = [ctd.Cruise ; repmat({'SKQ201818S'},4,1)];
    ctd.Cast     = [ctd.Cast;    repmat(2,4,1)];
    ctd.dnum     = [ctd.dnum;    repmat(737278.26318287,4,1)];
    ctd.Station  = [ctd.Station; repmat({'CEO_CTD'},4,1)];
    ctd.density  = [ctd.density; d'];
    ctd.pressure = [ctd.pressure;p'];
    ctd.salinity = [ctd.salinity;s'];
  case  {'CEO_2018' 'CEO_2018_2019'}
    % HLY1801 [https://web.whoi.edu/healy-1801/ctd-survey/]
    ctd_file = fullfile(cfg.path.calcasts,'2018','HLY1801_048.dcc');
    ctd1 = readtable(ctd_file,'FileType','text','HeaderLines',2,'ReadVariableNames',1);
    ctd1.Cast = repmat(48,size(ctd1.Pres));
    ctd1.Station = cellstr(repmat('UAF-CEO',size(ctd1.Pres)));
    ctd1.dnum    = repmat(datenum(2018,08,15,21,27,00),size(ctd1.Pres));
    
    ctd_file2 = fullfile(cfg.path.calcasts,'2018','HLY1801_049.dcc');
    ctd2 = readtable(ctd_file2,'FileType','text','HeaderLines',2,'ReadVariableNames',1);
    ctd2.Cast = repmat(49,size(ctd2.Pres));
    ctd2.Station = cellstr(repmat('DBO4-5N',size(ctd2.Pres)));
    ctd2.dnum    = repmat(datenum(2018,08,15,22,05,00),size(ctd2.Pres));
   
    % merge both casts into single table
    ctd_hly = [ctd1; ctd2];
    % calculate density, datenum, and ctd1
    ctd_hly.density  = sw_dens(ctd_hly.Sal_1_,ctd_hly.T90_1_,ctd_hly.Pres);
    ctd_hly.pressure = ctd_hly.Pres;
    ctd_hly.Cruise   = cellstr(repmat('HLY1801',size(ctd_hly.pressure)));
    ctd_hly.salinity = ctd_hly.Sal_1_;
    
    % Add OS1901 Cast: #33 Station: CEO
    file1 = fullfile(cfg.path.calcasts,'OS1901','os1901l1c033_ctd.nc');
    ctd.pressure = squeeze(ncread(file1,'P_1'));     % [DB]
    ctd.salinity = squeeze(ncread(file1,'S_41'));    % [PSU]
    ctd.density = sw_dens(ctd.salinity,squeeze(ncread(file1,'T_28')),ctd.pressure); % calculate density with seawater toolbox
    ctd.Cruise  = cellstr(repmat('OS2019',size(ctd.density)));
    ctd.Cast    = repmat(33,size(ctd.density));
    ctd.dnum    = repmat(datenum(2019,08,19,14,14,00),size(ctd.density));
    ctd.Station = cellstr(repmat('CEO',size(ctd.density)));
    
    % Combine all casts
    ctd.pressure = [ctd.pressure; ctd_hly.pressure];
    ctd.salinity = [ctd.salinity; ctd_hly.salinity];
    ctd.Cruise   = [ctd.Cruise; ctd_hly.Cruise];
    ctd.Cast     = [ctd.Cast; ctd_hly.Cast];
    ctd.dnum     = [ctd.dnum; ctd_hly.dnum];
    ctd.Station  = [ctd.Station; ctd_hly.Station];
    ctd.density  = [ctd.density; ctd_hly.density];
    
  case 'CEO_2019'
    % OS1901 Cast 33 | CEO
    file1 = fullfile(cfg.path.calcasts,'OS1901','os1901l1c033_ctd.nc');
    ctd.pressure = squeeze(ncread(file1,'P_1'));     % [DB]
    ctd.salinity = squeeze(ncread(file1,'S_41'));    % [PSU]
    ctd.density = sw_dens(ctd.salinity,squeeze(ncread(file1,'T_28')),ctd.pressure); % calculate density with seawater toolbox
    ctd.Cruise  = cellstr(repmat('OS2019',size(ctd.density)));
    ctd.Cast    = repmat(33,size(ctd.density));
    ctd.dnum    = repmat(datenum(2019,08,19,14,14,00),size(ctd.density));
    ctd.Station = cellstr(repmat('CEO',size(ctd.density)));
    
    % % No SUNA was deployed on the Norseman II 2020 cruise... 
    %     % Noraseman II 2020 Cast 1 | CEO
    %     file2 = fullfile(cfg.path.calcasts,'NS2020','DBO_NS20_ctd_L2_v1.csv');
    %     T = readtable(file2);
    %     T = T(strcmp(T.Station,'CEO'),:);
    %     ctd.pressure = [ctd.pressure; T.Pressure__dbar_];
    %     ctd.salinity = [ctd.salinity; T.Salinity__psu_];
    %     ctd.Cruise   = [ctd.Cruise; T.Cruise];
    %     ctd.Cast     = [ctd.Cast; T.Cast];
    %     ctd.dnum     = [ctd.dnum; datenum(T.Date_Time,'yyyy-mm-ddTHH:MM:SS')];
    %     ctd.Station  = [ctd.Station; T.Station];
    %     T.density = sw_dens(T.Salinity__psu_,T.Temperature__C_,T.Pressure__dbar_);% calculate density with seawater toolbox
    %     ctd.density  = [ctd.density; T.density];
    
  otherwise
    fprintf('need to set up read full ctd profile(s)\n')
    keyboard
end

%% Add unique identifier by combining cruise and station name
% In case cruises have the same station name, combine them into a unique
% identifier.
ctd.cruise_station = strcat(ctd.Cruise,'_',ctd.Station);
ctd.cruise_station = strrep(ctd.cruise_station,'-','_');
ctd.cruise_station = strrep(ctd.cruise_station,' ','');
ctd.cruise_station = strrep(ctd.cruise_station,'.','_');
ctd.cruise_station = strrep(ctd.cruise_station,'/','_');


end %% MAIN FUNCTION