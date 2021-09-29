function [cfg,offset] = SUNA_calibration_discrete_samples_all
%FUNCTION SUNA_calibration_discrete_samples
%
%  Syntax:
%    [data_proc,meta_proc,cfg] = SUNA_calibration_discrete_samples(data_proc,meta_proc,cfg)
%
%  Description: 
%
%  Syntax:
%
%  Description:
%
%  Examples:
%
%  References:
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Initialize script variables
dbstop if error
%% Set up script
mooring_lat = 71.6;
mooring_lon = -161.52;
% Read all calibration cast data
cfg.project = ['CEO_2015'];
cfg.year = 2015;
cfg.yeardir = ['F:\CEO\DataByDeploymentYear\2015\'];
cfg.datadir = fullfile(cfg.yeardir,'SUNA',filesep);
cfg = CEO_config(cfg);
cfg = CEO_read_calibration_cast_data(cfg,0); % the 0 here indicates do NOT limit to single cruise
%% Load SUNA data
% Automated quality assurance and quality control tests following best
% practices (where available)
% if ~exist('d2015','var')
%   d2015 = load('F:\CEO\DataByDeploymentYear\2015\SUNA\CEO_2015_SUNA_proc.mat','segments');
%   d2016 = load('F:\CEO\DataByDeploymentYear\2016\SUNA\CEO_2016_SUNA_proc.mat','segments');
%   d2017 = load('F:\CEO\DataByDeploymentYear\2017\SUNA\CEO_2017_SUNA_proc.mat','segments');
%   d2018 = load('F:\CEO\DataByDeploymentYear\2018\SUNA\CEO_2018_SUNA_proc.mat','segments');
%   d2017.segments.project = repmat({'CEO_2017'},d2017.segments.number_of_segments,1);
%   d2016.segments.project = repmat({'CEO_2016'},d2016.segments.number_of_segments,1);
%   d2015.segments.project = repmat({'CEO_2015'},d2015.segments.number_of_segments,1);
%   d2018.segments.project = repmat({'CEO_2018'},d2018.segments.number_of_segments,1);
%   data.project  = [d2015.segments.project; d2016.segments.project; d2017.segments.project; d2018.segments.project];
%   data.datenum  = [d2015.segments.datenum; d2016.segments.datenum; d2017.segments.datenum; d2018.segments.datenum];
%   data.Pressure = [d2015.segments.pressure; d2016.segments.pressure; d2017.segments.pressure; d2018.segments.pressure];
%   data.density  = [d2015.segments.density; d2016.segments.density; d2017.segments.density; d2018.segments.density];
%   data.Nitrate_corrected      = [d2015.segments.NO3_uM_cor_TCSS; d2016.segments.NO3_uM_cor_TCSS; d2017.segments.NO3_uM_cor_TCSS; d2018.segments.NO3_uM_cor_TCSS];
%   data.flag     = [d2015.segments.flag; d2016.segments.flag; d2017.segments.flag; d2018.segments.flag];
% end

for nyear = 2015:2018
  file = fullfile('F:\CEO\DataByDeploymentYear\',num2str(nyear),'SUNA',['CEO2_' num2str(nyear) '-' num2str(nyear+1) '_SUNAv2_Nitrate_L3_v1.csv']);
  opt = detectImportOptions(file);
  opt.VariableUnitsLine = opt.VariableNamesLine+1;
  
  if nyear == 2015
    data = readtable(file,opt);
    data.project = repmat({'CEO_2015'},size(data.Nitrate_raw,1),1);
  else
     dat = readtable(file,opt);
     dat.project = repmat({['CEO_' num2str(nyear)]},size(dat.Nitrate_raw,1),1);
     try
       data = [data; dat];
     catch
       fprintf('stopped here...\n');
       keyboard
     end
  end
end

% Calculate density
d.p = data.Pressure;
d.t = data.Temperature;
d.Sp = data.Salinity;
d.lat = mooring_lat;
d.lon = mooring_lon;
dens = gsw_rho_irving(d);
data.density = dens.rho;
data.datetime = datetime(strrep(data.DateTime,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:SS');
data.datenum  = datenum(data.datetime);
%% Calculate distance in time and space (horizontally and vertically)
cfg.calcasts.distkm = nan(size(cfg.calcasts.CastID));
cfg.calcasts.distdb = nan(size(cfg.calcasts.CastID));
cfg.calcasts.disthr = nan(size(cfg.calcasts.CastID));
for nsamp = 1:size(cfg.calcasts,1)
  % find distance horizontally and vertically
  cfg.calcasts.distkm(nsamp) = round(m_lldist([cfg.calcasts.Lon(nsamp) mooring_lon],[cfg.calcasts.Lat(nsamp) mooring_lat]),2);     % km separating measurements
  % Calculate pressure
  if isnan(cfg.calcasts.Pressure(nsamp))
    cfg.calcasts.Pressure(nsamp) = gsw_p_from_z(-cfg.calcasts.DepSM(nsamp),cfg.calcasts.Lat(nsamp));
  end
  cfg.calcasts.distdb(nsamp) = round(abs(cfg.calcasts.Pressure(nsamp)-nanmean(data.Pressure)),2);
  % find distance in time
  time_dist = round(cfg.calcasts.dnum(nsamp)-data.datenum,2);
  idx_time_near = nearest(time_dist,0); % distance in days
  cfg.calcasts.disthr(nsamp) = round(abs(etime(datevec(data.datenum(idx_time_near)),datevec(cfg.calcasts.dnum(nsamp))))./60./60,2);
  %deploy_hours  = round(abs(etime(datevec(cfg.mooring.deploydate),datevec(cfg.discrete_ref.dnum(nsamp))))./60./60,2);  % Hours separating measurements (2 decimal places)
  %recover_hours = round(abs(etime(datevec(cfg.mooring.recoverdate),datevec(cfg.discrete_ref.dnum(nsamp))))./60./60,2); % Hours separating measurements from recovery
  % Calculate density
  if isnan(cfg.calcasts.density(nsamp))
    cfg.calcasts.density(nsamp) = sw_dens(cfg.calcasts.Salinity(nsamp),cfg.calcasts.Temp(nsamp), cfg.calcasts.Pressure(nsamp));
  end
end

%% Calibration casts
cfg.calcasts.NO3_variable = cell(size(cfg.calcasts.CastID));
% Loop through castID
u_casts = unique(cfg.calcasts.CastID,'stable');
idx_match = logical(zeros(size(cfg.calcasts.CastID)));

for nsamp = 1:size(cfg.calcasts,1)
  cal = cfg.calcasts(nsamp,:);
  fprintf('Looking for nearby discrete samples | %s %s\n',char(cal.CRUISE), char(cal.StationID))
  % match micromolar NO3
  if isfinite(cal.NO3_uM)
    idx_match(nsamp) = true;
    cal.NO3_variable = {'NO3_uM'};
  elseif isfinite(cal.NO3) && ~isfinite(cal.NO3_uM)
    fprintf('  NO3 umol/kg values found... converting to umol/L!\n')
    idx_match(nsamp) = true;
    if isnan(cal.Pressure)
      cal.Pressure = gsw_p_from_z(-1.0*cal.DepSM,cal.Lat);
      dens = sw_dens(cal.Salinity,cal.Temp, cal.Pressure);
    else
       dens = sw_dens(cal.Salinity,cal.Temp, cal.Pressure);
    end
    cal.NO3_uM = round(cal.NO3 .* 1000 ./dens,3); % [umol/L] * 1000[L/m3] * 1/[kg/m3]
    cal.NO3_variable = {'NO3_uM'};
  elseif isfinite(cal.NO3_NO2_uM)
    fprintf('  NO3 + NO2 uM values found\n')
    idx_match(nsamp) = true;
    cal.NO3_variable = {'NO3_NO2_uM'};
  else
    %fprintf('No Nitrate or nitrite values found...\n')
  end
  try
    cfg.calcasts(nsamp,:) = cal;
  catch
    fprintf('stopped here\n')
    keyboard
  end
end
calcasts = cfg.calcasts(idx_match == 1,:);
calcasts(calcasts.distdb > 7,:) = [];


%% PLOT
years = 2015:2018;
clrs = lansey(numel(years)); % brewermap(numel(years),'RdYlBu'); %
makefig; ax = gca;  hold(ax,'on'); grid(ax,'on'); ax.YDir = 'normal'; ax.Color = [0.75 0.75 .75];
u_project = unique(data.project);
for nn = 1:numel(years)
  idx_project = strcmp(data.project,u_project{nn});
  plot(ax,data.datenum(idx_project),data.Nitrate_corrected(idx_project),'o','Color',clrs(nn,:),'MarkerSize',7,'MarkerFaceColor',clrs(nn,:),'DisplayName',strrep(u_project{nn},'_','\_'))
end

clrs = jet(numel(calcasts.dnum));
for nd = 1:numel(calcasts.dnum)
  no3_var = calcasts.NO3_variable{nd};
  plot_text = [pad(calcasts.CRUISE{nd},10) ' | Station: ' pad(calcasts.StationID{nd},10) ' Cast#'  num2str(calcasts.CastID(nd)) ...
    ' | ' strrep(no3_var,'_','\_') ...
    ' @ ' num2str(calcasts.distdb(nd),3) 'dbar off' ...
    ' | ' num2str(calcasts.distkm(nd)) 'km away' ...
    ' & ' num2str(calcasts.disthr(nd)) 'hr off'];
  plot(ax,calcasts.dnum(nd),calcasts.(no3_var)(nd),'kh','MarkerSize',15,'LineWidth',0.2,'MarkerFaceColor',clrs(nd,:),'DisplayName',plot_text);
end

hl = legend(ax,'show'); hl.FontSize = 10; hl.Location = 'northwest';
datetick(ax,'x');
axis(ax,'tight');

%% 
makefig; ax = gca;  hold(ax,'on'); grid(ax,'on'); ax.YDir = 'rev'; ax.Color = [0.75 0.75 .75];
scatter(ax,data.datenum,data.Pressure,30,data.Nitrate_corrected,'filled','displayName','Nitrate corrected');
gd = isfinite(cfg.calcasts.NO3_uM);
scatter(ax,cfg.calcasts.dnum(gd),cfg.calcasts.Pressure(gd),100,cfg.calcasts.NO3_uM(gd),'h','filled','MarkerEdgeColor','k','DisplayName','Discrete NO_3 uM')
gd = isfinite(cfg.calcasts.NO3_NO2_uM);
scatter(ax,cfg.calcasts.dnum(gd),cfg.calcasts.Pressure(gd),100,cfg.calcasts.NO3_NO2_uM(gd),'h','filled','MarkerEdgeColor','k','DisplayName','Discrete NO_3+NO_2 uM')
hl = legend(ax,'show'); hl.FontSize = 10; 
datetick(ax,'x');
ax.YLim = [25 45];
cb = colorbar(ax);
cb.Label.String = 'Nitrate [\muM]';
ax.YLabel.String = 'Pressure [dbar]';
cb.Limits = [0 ceil(nanmax(data.Nitrate_corrected))];
ax.CLim   = cb.Limits;
% Plot interpolated values
p = 1:45;
for nn = 1:numel(u_casts)
  idx  = cfg.calcasts.CastID == u_casts(nn);
  cast = cfg.calcasts(idx,:);
  uname= unique(cast.StationID,'stable');
  for ii = 1:numel(uname)
    ncast = cast(strcmp(cast.StationID,uname{ii}),:);
    no3_var = ncast.NO3_variable{1};
    if ~isempty(no3_var)
      plot_text = [pad(ncast.CRUISE{1},10) ' | Station: ' pad(ncast.StationID{1},10) ' Cast#'  num2str(ncast.CastID(1)) ...
        ' | ' strrep(no3_var,'_','\_') ...
        ' | ' num2str(ncast.distkm(1)) 'km away' ...
        ' & ' num2str(ncast.disthr(1)) 'hr off'];
      
      no3_interp = interp1(ncast.Pressure,ncast.(no3_var),p);
      dnum = repmat(ncast.dnum(1),size(no3_interp));
      gd = isfinite(no3_interp);
      scatter(ax,dnum(gd),p(gd),30,no3_interp(gd),'>','filled','MarkerEdgeColor','k','DisplayName',['Interpolated: ' plot_text])
    end
  end
end

keyboard
end %% MAIN FUNCTION

