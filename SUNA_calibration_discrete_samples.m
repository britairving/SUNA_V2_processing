function [cfg,offset] = SUNA_calibration_discrete_samples(data_avg,meta_proc,cfg)
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
%% Set up script
if isfield(cfg,'save_figures')
  save_fig = cfg.save_figures;
else
  save_fig = false;
end
close all;

%% Calculate distance in time and space (horizontally and vertically)
cfg.calcasts.distkm = nan(size(cfg.calcasts.CastID));
cfg.calcasts.distdb = nan(size(cfg.calcasts.CastID));
cfg.calcasts.disthr = nan(size(cfg.calcasts.CastID));
for nsamp = 1:size(cfg.calcasts,1)
  % find distance horizontally and vertically
  cfg.calcasts.distkm(nsamp) = round(m_lldist([cfg.calcasts.Lon(nsamp) cfg.mooring.longitude],[cfg.calcasts.Lat(nsamp) cfg.mooring.latitude]),2);     % km separating measurements
  % Calculate pressure
  if isnan(cfg.calcasts.Pressure(nsamp))
    cfg.calcasts.Pressure(nsamp) = gsw_p_from_z(-cfg.calcasts.DepSM(nsamp),cfg.calcasts.Lat(nsamp));
  end
  cfg.calcasts.distdb(nsamp) = round(abs(cfg.calcasts.Pressure(nsamp)-nanmean(data_avg.pressure)),2);
  % find distance in time
  time_dist = round(cfg.calcasts.dnum(nsamp)-data_avg.datenum,2);
  idx_time_near = nearest(time_dist,0); % distance in days
  cfg.calcasts.disthr(nsamp) = round(abs(etime(datevec(data_avg.datenum(idx_time_near)),datevec(cfg.calcasts.dnum(nsamp))))./60./60,2);
  %deploy_hours  = round(abs(etime(datevec(cfg.mooring.deploydate),datevec(cfg.discrete_ref.dnum(nsamp))))./60./60,2);  % Hours separating measurements (2 decimal places)
  %recover_hours = round(abs(etime(datevec(cfg.mooring.recoverdate),datevec(cfg.discrete_ref.dnum(nsamp))))./60./60,2); % Hours separating measurements from recovery
  % Calculate density
  if isnan(cfg.calcasts.density(nsamp))
    cfg.calcasts.density(nsamp) = sw_dens(cfg.calcasts.Salinity(nsamp),cfg.calcasts.Temp(nsamp), cfg.calcasts.Pressure(nsamp));
  end
end

%% Calibration casts
cfg.calcasts.NO3_variable = cell(size(cfg.calcasts.CastID));
cfg.calcasts.cal_match    = logical(zeros(size(cfg.calcasts.CastID)));
% Loop through castID
u_casts = unique(cfg.calcasts.CastID,'stable');
idx_match = logical(zeros(size(cfg.calcasts.CastID)));
for ncast = 1:numel(u_casts)
  idx_cast = find(cfg.calcasts.CastID == u_casts(ncast));
  
  fprintf('Looking for nearby discrete samples | %s %s\n',cfg.calcasts.CRUISE{idx_cast(1)}, cfg.calcasts.StationID{idx_cast(1)})
  % match micromolar NO3
  if any(isfinite(cfg.calcasts.NO3_uM(idx_cast)))
    idx_match(idx_cast) = true;
    cfg.calcasts.NO3_variable(idx_cast) = {'NO3_uM'};
    cfg.calcasts.cal_match(idx_cast) = 1;
  elseif any(isfinite(cfg.calcasts.NO3(idx_cast))) && ~any(isfinite(cfg.calcasts.NO3_uM(idx_cast)))
    fprintf('  NO3 umol/kg values found... converting to umol/L!\n')
    %if all(isfinite( cfg.calcasts.Sigmat11(idx_cast)))
    %  dens = (cfg.calcasts.Sigmat11(idx_cast) + 1000.0); % add 1000 to SigmaT to get density kg/m3
    %else
    if all(isnan(cfg.calcasts.Pressure(idx_cast)))
      cfg.calcasts.Pressure(idx_cast) = gsw_p_from_z(-1.0*cfg.calcasts.DepSM(idx_cast),cfg.calcasts.Lat(idx_cast));
      dens = sw_dens(cfg.calcasts.Salinity(idx_cast),cfg.calcasts.Temp(idx_cast), cfg.calcasts.Pressure(idx_cast));
    else
      dens = sw_dens(cfg.calcasts.Salinity(idx_cast),cfg.calcasts.Temp(idx_cast), cfg.calcasts.Pressure(idx_cast));
    end
    %end
    cfg.calcasts.NO3_uM(idx_cast) = round(cfg.calcasts.NO3(idx_cast) .* 1000 ./dens,3); % [umol/L] * 1000[L/m3] * 1/[kg/m3]
    idx_match(idx_cast) = true;
    cfg.calcasts.NO3_variable(idx_cast) = {'NO3_uM'};
    cfg.calcasts.cal_match(idx_cast) = 1;
  elseif any(isfinite(cfg.calcasts.NO3_NO2_uM(idx_cast)))
    fprintf('  NO3 + NO2 uM values found\n')
    idx_match(idx_cast) = true;
    cfg.calcasts.NO3_variable(idx_cast) = {'NO3_NO2_uM'};
    cfg.calcasts.cal_match(idx_cast) = 1;
  else
    %fprintf('No Nitrate or nitrite values found...\n')
  end
end
[casts,iu] = unique(cfg.calcasts.StationID(cfg.calcasts.cal_match),'stable');


%% Read full depth ctd profile to see if okay to use...
switch cfg.project
  case 'CEO_2016'
    % AMBON2017
    ctd_file = fullfile(cfg.path.calcasts,'AMBON2015_2017','AMBON2017_ctd_L3_v1.csv');
    ctd = readtable(ctd_file);
    ctd_match = ismember(ctd.Station,unique(cfg.calcasts.StationID));  %ctd_match = strcmp(ctd.Station,cfg.discrete_ref.StationID{1});
    ctd = ctd(ctd_match,:);
    % calculate density, datenum, and datetime
    ctd.density = sw_dens(ctd.salinity__psu_,ctd.temperature__C_,ctd.pressure__dbar_);
    ctd.dtime = datetime(strrep(ctd.Date_Time,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:SS');
    ctd.dnum  = datenum(ctd.dtime);
    ctd.pressure = ctd.pressure__dbar_;
    ctd.salinity = ctd.salinity__psu_;
  case 'CEO_2015'
    % AMBON2015
    ctd_file = fullfile(cfg.path.calcasts,'AMBON2015_2017','AMBON2015_ctd_L3_v1.csv');
    ctd = readtable(ctd_file);
    ctd_match = ismember(ctd.Station,unique(cfg.calcasts.StationID));  %ctd_match = strcmp(ctd.Station,cfg.discrete_ref.StationID{1});
    ctd = ctd(ctd_match,:);
    % calculate density, datenum, and datetime
    ctd.density = sw_dens(ctd.salinity__psu_,ctd.temperature__C_,ctd.pressure__dbar_);
    ctd.dtime = datetime(strrep(ctd.Date_Time,'T',' '),'InputFormat','yyyy-MM-dd HH:mm:SS');
    ctd.dnum  = datenum(ctd.dtime);
    ctd.pressure = ctd.pressure__dbar_;
    ctd.salinity = ctd.salinity__psu_;
  case 'CEO_2017'
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
    ctd_match = ismember(ctd3.Station,unique(cfg.calcasts.StationID));  %ctd_match = strcmp(ctd3.Station,cfg.discrete_ref.StationID{1});
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
  case  'CEO_2018'
    % HLY1801 [https://web.whoi.edu/healy-1801/ctd-survey/]
    ctd_file = fullfile(cfg.path.calcasts,'2018','HLY1801_048.dcc');
    ctd = readtable(ctd_file,'FileType','text','HeaderLines',2,'ReadVariableNames',1);
    ctd.Cast = repmat(48,size(ctd.Pres));
    ctd.Station = cellstr(repmat('UAF-CEO',size(ctd.Pres)));
    ctd.dnum    = repmat(datenum(2018,08,15,21,27,00),size(ctd.Pres));
    
    ctd_file2 = fullfile(cfg.path.calcasts,'2018','HLY1801_049.dcc');
    ctd2 = readtable(ctd_file2,'FileType','text','HeaderLines',2,'ReadVariableNames',1);
    ctd2.Cast = repmat(49,size(ctd2.Pres));
    ctd2.Station = cellstr(repmat('DBO4-5N',size(ctd2.Pres)));
    ctd2.dnum    = repmat(datenum(2018,08,15,22,05,00),size(ctd2.Pres));
    % merge both casts into single table
    ctd = [ctd; ctd2];
    % calculate density, datenum, and ctd
    ctd.density  = sw_dens(ctd.Sal_1_,ctd.T90_1_,ctd.Pres);
    ctd.pressure = ctd.Pres;
    ctd.Cruise   = cellstr(repmat('HLY1801',size(ctd.pressure)));
    ctd.salinity = ctd.Sal_1_;
  otherwise
    fprintf('need to set up read full ctd profile(s)\n')
    keyboard
end

%% pull out SUNA and CTD Data in closest casts
mtch = struct();
for ncast = 1:numel(casts)
  cast_name = strrep(casts{ncast},'.','_');
  cast_name = strrep(cast_name,'-','_');
  mtch.(cast_name) = struct();
  mtch.(cast_name).idx_ctd = find(strcmp(ctd.Station,casts{ncast}));
  mtch.(cast_name).idx_bot = find(strcmp(cfg.calcasts.StationID,casts{ncast}));
  % only include those with good data
  no3_var = cfg.calcasts.NO3_variable{ mtch.(cast_name).idx_bot(1)};
  mtch.(cast_name).idx_bot = mtch.(cast_name).idx_bot(isfinite(cfg.calcasts.(no3_var)(mtch.(cast_name).idx_bot)));
end
% now find 2 day range for SUNA data
for ncast = 1:numel(casts)
  cast_name = strrep(casts{ncast},'.','_');
  cast_name = strrep(cast_name,'-','_');
  trng = find(data_avg.datenum >= ctd.dnum(mtch.(cast_name).idx_ctd(1))-2 ...
    & data_avg.datenum <= ctd.dnum(mtch.(cast_name).idx_ctd(1))+2 ...
    & data_avg.flag    <= meta_proc.flag.not_evaluated);
  % now fill mtch struture
  mtch.(cast_name).idx_SUNA         = trng;
  mtch.(cast_name).SUNA_pressure    = round(nanmean(data_avg.pressure(trng)),2);
  mtch.(cast_name).SUNA_salinity    = round(nanmean(data_avg.salinity(trng)),2);
  mtch.(cast_name).SUNA_density     = round(nanmean(data_avg.density(trng)),2);
  mtch.(cast_name).SUNA_NO3_uM      = round(nanmean(data_avg.NO3_uM_smo(trng)),2);
  mtch.(cast_name).SUNA_NO3_uM_TCSS = round(nanmean(data_avg.NO3_uM_TCSS_smo(trng)),2);
  if isfield(data_avg,'NO3_uM_TCSS_cor')
    mtch.(cast_name).SUNA_NO3_uM_TCSS_cor = round(nanmean(data_avg.NO3_uM_TCSS_cor(trng)),2);
  end
end

%% Plot1 | CTD cast with locations of discrete samples
makefig; ax = subplot(1,5,1:3); ax2 = subplot(1,5,4); ax3 = subplot(1,5,5);
hold(ax,'on');  grid(ax,'on');  set(ax,'YDir','rev');
hold(ax2,'on'); grid(ax2,'on'); set(ax2,'YDir','rev');
hold(ax3,'on'); grid(ax3,'on'); set(ax3,'YDir','rev');
ax.XLabel.String = 'Density [kg/m^3]';
ax.YLabel.String = 'Pressure [dbar]';
ax2.XLabel.String = 'NO_3 [\muM]';
ax2.YLabel.String = 'Pressure [dbar]';
ax3.XLabel.String = 'Salinity';
ax3.YLabel.String = 'Pressure [dbar]';


clrs = jet(numel(casts));
for ncast = 1:numel(casts)
  cast_name = strrep(casts{ncast},'.','_');
  cast_name = strrep(cast_name,'-','_');
  idx_ctd = mtch.(cast_name).idx_ctd;
  idx_bot = mtch.(cast_name).idx_bot;
  ctd_str = [ctd.Station{idx_ctd(1)} ': Cast#' num2str(ctd.Cast(idx_ctd(1)))];
  bot_str = [ctd_str ' ' strrep(cfg.calcasts.NO3_variable{idx_bot(1)},'_','\_') ];
  no3_var = cfg.calcasts.NO3_variable{idx_bot(1)};
  imin = nearest(cfg.calcasts.Pressure(idx_bot),mtch.(cast_name).SUNA_pressure);
  imin = idx_bot(imin);
  plot_text = [pad(cfg.calcasts.CRUISE{imin},10) ' | Station: ' pad(cfg.calcasts.StationID{imin},10) ' Cast#'  num2str(cfg.calcasts.CastID(imin)) ...
    ' | ' strrep(no3_var,'_','\_') ...
    ...%' @ ' num2str(cfg.calcasts.Pressure(imin),3) ' dbar' ...
    ' | ' num2str(cfg.calcasts.distkm(imin)) 'km away' ...
    ' & ' num2str(cfg.calcasts.disthr(imin)) 'hr off'];
  %% Density vs Pressure
  % Plot full cast density
  plot(ax,ctd.density(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',plot_text);
  % Plot discrete bottle density
  plot(ax,cfg.calcasts.density(idx_bot),cfg.calcasts.Pressure(idx_bot),'kh','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' ctd_str]);
  % plot SUNA density
  plot(ax,mtch.(cast_name).SUNA_density,mtch.(cast_name).SUNA_pressure,'kd','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['SUNA Density @' ctd_str]);
  %% Nitrate vs Pressure
  plot(ax2,cfg.calcasts.(no3_var)(idx_bot),cfg.calcasts.Pressure(idx_bot),'k-h','MarkerFaceColor',clrs(ncast,:),'MarkerSize',10,'DisplayName',bot_str);
  plot(ax2,mtch.(cast_name).SUNA_NO3_uM,      mtch.(cast_name).SUNA_pressure,'ko','MarkerSize',10,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['SUNA NO_{3 raw} @' ctd_str]);
  plot(ax2,mtch.(cast_name).SUNA_NO3_uM_TCSS, mtch.(cast_name).SUNA_pressure,'kd','MarkerSize',10,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['SUNA NO_{3 TCSS} @ ' ctd_str]);
  if isfield(data_avg,'NO3_uM_TCSS_cor')
    plot(ax2,mtch.(cast_name).SUNA_NO3_uM_TCSS_cor, mtch.(cast_name).SUNA_pressure,'k>','MarkerSize',10,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['SUNA NO_{3 TCSS,COR} @ ' ctd_str]);
  end
  %% Salinity vs Pressure
  % Plot full cast density
  plot(ax3,ctd.salinity(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor',clrs(ncast,:),'MarkerSize',8,'DisplayName',plot_text);
  % Plot discrete bottle density
  plot(ax3,cfg.calcasts.Salinity(idx_bot),cfg.calcasts.Pressure(idx_bot),'kh','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' ctd_str]);
  % plot SUNA density
  plot(ax3,mtch.(cast_name).SUNA_salinity,mtch.(cast_name).SUNA_pressure,'kd','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['SUNA Salinity @' ctd_str]);
  
end


% Plot AMBON CEO16 and CEO17 profiles.... no nitrate but just for
% comparison sake
if strcmp(cfg.project,'CEO_2016')
  
  idx_bot = find(strcmp(cfg.calcasts.StationID,'CEO16'),1);
  idx_ctd = find(strcmp(ctd.Station,'CEO16'));
  plot_text = [pad(cfg.calcasts.CRUISE{idx_bot},10) ' | Station: ' pad(cfg.calcasts.StationID{idx_bot},10) ' Cast#'  num2str(cfg.calcasts.CastID(idx_bot)) ...
    ' | ' num2str(cfg.calcasts.distkm(idx_bot)) 'km away' ...
    ' & ' num2str(cfg.calcasts.disthr(idx_bot)) 'hr off'];
  plot(ax,ctd.density(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor','k','MarkerSize',8,'DisplayName',plot_text);
  
  idx_bot = find(strcmp(cfg.calcasts.StationID,'CEO17'),1);
  idx_ctd = find(strcmp(ctd.Station,'CEO17'));
  plot_text = [pad(cfg.calcasts.CRUISE{idx_bot},10) ' | Station: ' pad(cfg.calcasts.StationID{idx_bot},10) ' Cast#'  num2str(cfg.calcasts.CastID(idx_bot)) ...
    ' | ' num2str(cfg.calcasts.distkm(idx_bot)) 'km away' ...
    ' & ' num2str(cfg.calcasts.disthr(idx_bot)) 'hr off'];
  plot(ax,ctd.density(idx_ctd),ctd.pressure(idx_ctd),'k-s','MarkerFaceColor','g','MarkerSize',8,'DisplayName',plot_text);
end

axis(ax2,'tight');
ax2.YLim = ax.YLim;
ax2.Box = 'on';
ax.Title.String  = [strrep(cfg.project,'_','\_') ' | Nearby CTD casts and discrete samples'];
hl = legend(ax,'show','Location','southwest');
hl.FontSize =  12;
hl = legend(ax2,'show','Location','best');
hl.FontSize =  12;
hl.Position(1:2) = [0.53 0.69];% 0.12 0.33];

if save_fig
  savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison1']);
  standard_printfig_highrespng(savename);
end

%% Plot2 | SUNA datenum vs Nitrate with discrete samples
idx_good = data_avg.flag <= meta_proc.flag.not_evaluated;
% basic plot to show corrected values
makefig; a = gca; a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
plot(a,data_avg.datenum(idx_good),data_avg.NO3_uM_smo(idx_good),'ko','MarkerSize',5,'MarkerFaceColor','k',        'DisplayName','SUNA NO_{3 Raw}')
plot(a,data_avg.datenum(idx_good),data_avg.NO3_uM_TCSS_smo(idx_good),'kd','MarkerSize',5,'MarkerFaceColor','m','DisplayName','SUNA NO_{3 TCSS}')
if isfield(data_avg,'NO3_uM_TCSS_cor')
  plot(a,data_avg.datenum(idx_good),data_avg.NO3_uM_TCSS_cor(idx_good),'k>','MarkerSize',5,'MarkerFaceColor','g','DisplayName','SUNA NO_{3 TCSS,COR}')
end
% Plot discrete samples
for nclosest = 1:numel(casts)
  cast_name = strrep(casts{nclosest},'.','_');
  cast_name = strrep(cast_name,'-','_');
  imin = nearest(cfg.calcasts.Pressure(mtch.(cast_name).idx_bot),mtch.(cast_name).SUNA_pressure);
  imin = mtch.(cast_name).idx_bot(imin);
  no3_var = cfg.calcasts.NO3_variable{imin};
  plot_text = [pad(cfg.calcasts.CRUISE{imin},12) ' | Station: ' pad(cfg.calcasts.StationID{imin},10) ' Cast#'  num2str(cfg.calcasts.CastID(imin)) ...
    ' | ' strrep(no3_var,'_','\_') ...
    ...%' @ ' num2str(cfg.calcasts.Pressure(imin),3) ' dbar' ...
    ' | ' num2str(cfg.calcasts.distkm(imin)) 'km away' ...
    ' & ' num2str(cfg.calcasts.disthr(imin)) 'hr off'];
  plot(a,cfg.calcasts.dnum(imin),cfg.calcasts.(no3_var)(imin),'kh','MarkerSize',16,'MarkerFaceColor',clrs(nclosest,:),'LineWidth',1,'DisplayName',plot_text);
end
ylabel(a,'NO_3 [\muM]'); ylim(a,[-10 60]);
axis(a,'tight'); datetick(a,'x','keeplimits');
a.Title.String  = [strrep(cfg.project,'_','\_') ' | SUNA ' cfg.mooring.SN.SUNA ];
hl = legend(a,'show','Location','northwest');
hl.FontSize = 13;

if save_fig
  savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison2']);
  standard_printfig_highrespng(savename);
end

%% Plot3 | Density vs Nitrate
% makefig; a = gca; a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
 makefig; a = subplot(2,1,1);a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
if strcmp(cfg.project,'CEO_2017')
  casts = casts(2:end);
  clrs  = clrs(2:end,:);
end
ypos = 0.96;
for ncast = 1:numel(casts)
  % cast indicies and legend labels
  cast_name = strrep(casts{ncast},'.','_');
  cast_name = strrep(cast_name,'-','_');
  idx_SUNA  = mtch.(cast_name).idx_SUNA;
  idx_bot = mtch.(cast_name).idx_bot;
  idx_bot = idx_bot(cfg.calcasts.Pressure(idx_bot) > 25);  %Only look at relationship in lower water column
  no3_var = cfg.calcasts.NO3_variable{idx_bot(1)};
  bot_str = [strrep(cast_name,'_','\_') ' ' strrep(no3_var,'_','\_')];
  
  % SUNA raw data
  plot(a,data_avg.density(idx_SUNA),data_avg.NO3_uM_smo(idx_SUNA),'ko','MarkerSize',10,'MarkerFaceColor','k', 'DisplayName','SUNA NO_{3 Raw}')
  [p,S] = polyfit(data_avg.density(idx_SUNA),data_avg.NO3_uM_smo(idx_SUNA),1); no3fit = polyval(p,data_avg.density(idx_SUNA),S);
  plot(a,data_avg.density(idx_SUNA),no3fit,'k-','LineWidth',2,'DisplayName','SUNA: Density NO_{3 Raw} regression')
  text(a,0.02,ypos,['NO_{3raw} Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  % SUNA TCSS corrected data
  plot(a,data_avg.density(idx_SUNA),data_avg.NO3_uM_TCSS_smo(idx_SUNA),'kd','MarkerSize',10,'MarkerFaceColor','m','DisplayName','SUNA NO_{3 TCSS}')
  [p,S] = polyfit(data_avg.density(idx_SUNA),data_avg.NO3_uM_TCSS_smo(idx_SUNA),1); no3fit = polyval(p,data_avg.density(idx_SUNA),S);
  plot(a,data_avg.density(idx_SUNA),no3fit,'m-','LineWidth',2,'DisplayName','SUNA: Density NO_{3 TCSS} regression')
  text(a,0.02,ypos-0.05,['NO_{3TCSS} Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  % SUNA TCSS corrected data AND corrected with relative offset from
  % calfile
  if isfield(data_avg,'NO3_uM_TCSS_cor')
    plot(a,data_avg.density(idx_SUNA),data_avg.NO3_uM_TCSS_cor(idx_SUNA),'k>','MarkerSize',10,'MarkerFaceColor','g','DisplayName','SUNA NO_{3 TCSS,COR}')
    [p,S] = polyfit(data_avg.density(idx_SUNA),data_avg.NO3_uM_TCSS_cor(idx_SUNA),1); no3fit = polyval(p,data_avg.density(idx_SUNA),S);
    plot(a,data_avg.density(idx_SUNA),no3fit,'g-','LineWidth',2,'DisplayName','SUNA: Density NO_{3 TCSS,COR} regression')
    text(a,0.02,ypos-0.1,['NO_{3TCSS,COR} Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  end
  
  % Discrete Samples
  plot(a,cfg.calcasts.density(idx_bot),cfg.calcasts.(no3_var)(idx_bot),'kh','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' bot_str]);
  % Discrete sample linear fit
  [p,S] = polyfit(cfg.calcasts.density(idx_bot),cfg.calcasts.(no3_var)(idx_bot),1);  no3fit = polyval(p,cfg.calcasts.density(idx_bot),S);
  plot(a,cfg.calcasts.density(idx_bot),no3fit,'--','MarkerSize',5,'Color',clrs(ncast,:),'LineWidth',2,'DisplayName','Discrete samples: Density NO_3 regression')
  text(a,0.02,ypos-0.15,[cast_name ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  ypos = ypos - numel(casts)*0.1;
end
a.XLabel.String = 'Density [kg/m^3]';
a.YLabel.String =  'NO_3 [\muM]';
a.Title.String  = [strrep(cfg.project,'_','\_') ' | Density - Nitrate linear regression'];
hl = legend(a,'show');
hl.FontSize = 12;
hl.Location = 'bestoutside';
% if save_fig
%   savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison3']);
%   standard_printfig_highrespng(savename);
% end

%% Plot4 | Salinity vs Nitrate
% makefig; a = gca; a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
a = subplot(2,1,2); a.Color = [0.75 0.75 0.75];hold(a,'on'); grid(a,'on'); a.YDir = 'normal';
ypos = 0.96;
for ncast = 1:numel(casts)
  % cast indicies and legend labels
  cast_name = strrep(casts{ncast},'.','_');
  cast_name = strrep(cast_name,'-','_');
  idx_SUNA  = mtch.(cast_name).idx_SUNA;
  idx_bot = mtch.(cast_name).idx_bot;
  idx_bot = idx_bot(cfg.calcasts.Pressure(idx_bot) > 25);  %Only look at relationship in lower water column
  no3_var = cfg.calcasts.NO3_variable{idx_bot(1)};
  bot_str = [strrep(cast_name,'_','\_') ' ' strrep(no3_var,'_','\_')];
  
  % SUNA raw data
  plot(a,data_avg.salinity(idx_SUNA),data_avg.NO3_uM_smo(idx_SUNA),'ko','MarkerSize',10,'MarkerFaceColor','k', 'DisplayName','SUNA NO_{3 Raw}')
  [p,S] = polyfit(data_avg.salinity(idx_SUNA),data_avg.NO3_uM_smo(idx_SUNA),1); no3fit = polyval(p,data_avg.salinity(idx_SUNA),S);
  plot(a,data_avg.salinity(idx_SUNA),no3fit,'k-','LineWidth',2,'DisplayName','SUNA: salinity NO_{3 Raw} regression')
  text(a,0.02,ypos,['NO_{3raw} Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  % SUNA TCSS corrected data
  plot(a,data_avg.salinity(idx_SUNA),data_avg.NO3_uM_TCSS_smo(idx_SUNA),'kd','MarkerSize',10,'MarkerFaceColor','m','DisplayName','SUNA NO_{3 TCSS}')
  [p,S] = polyfit(data_avg.salinity(idx_SUNA),data_avg.NO3_uM_TCSS_smo(idx_SUNA),1); no3fit = polyval(p,data_avg.salinity(idx_SUNA),S);
  plot(a,data_avg.salinity(idx_SUNA),no3fit,'m-','LineWidth',2,'DisplayName','SUNA: salinity NO_{3 TCSS} regression')
  text(a,0.02,ypos-0.05,['NO_{3TCSS} Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  % SUNA TCSS corrected data AND corrected with relative offset from
  % calfile
  if isfield(data_avg,'NO3_uM_TCSS_cor')
    plot(a,data_avg.salinity(idx_SUNA),data_avg.NO3_uM_TCSS_cor(idx_SUNA),'k>','MarkerSize',10,'MarkerFaceColor','g','DisplayName','SUNA NO_{3 TCSS,COR}')
    [p,S] = polyfit(data_avg.salinity(idx_SUNA),data_avg.NO3_uM_TCSS_cor(idx_SUNA),1); no3fit = polyval(p,data_avg.salinity(idx_SUNA),S);
    plot(a,data_avg.salinity(idx_SUNA),no3fit,'g-','LineWidth',2,'DisplayName','SUNA: salinity NO_{3 TCSS,COR} regression')
    text(a,0.02,ypos-0.1,['NO_{3TCSS,COR} Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  end
  
  % Discrete Samples
  plot(a,cfg.calcasts.Salinity(idx_bot),cfg.calcasts.(no3_var)(idx_bot),'kh','MarkerSize',12,'MarkerFaceColor',clrs(ncast,:),'LineWidth',1,'DisplayName',['Discrete samples ' bot_str]);
  % Discrete sample linear fit
  [p,S] = polyfit(cfg.calcasts.Salinity(idx_bot),cfg.calcasts.(no3_var)(idx_bot),1);  no3fit = polyval(p,cfg.calcasts.Salinity(idx_bot),S);
  plot(a,cfg.calcasts.Salinity(idx_bot),no3fit,'--','MarkerSize',5,'Color',clrs(ncast,:),'LineWidth',2,'DisplayName','Discrete samples: Salinity NO_3 regression')
  text(a,0.02,ypos-0.15,[cast_name ' Slope = ' num2str(p(1),'%.2f')],'units','normalized','FontSize',14,'FontWeight','bold','Color',clrs(ncast,:))
  ypos = ypos - numel(casts)*0.1;
end
a.XLabel.String = 'Salinity';
a.YLabel.String =  'NO_3 [\muM]';
a.Title.String  = 'Salinity - Nitrate linear regression';
hl = legend(a,'show');
hl.FontSize = 12;
hl.Location = 'bestoutside';
if save_fig
  savename = fullfile(cfg.datadir,[cfg.project '_discretebottle_comparison3']);
  standard_printfig_highrespng(savename);
end

%% FINALLY - CHOOSE WHAT TO DO
DONE_CHOOSING = 0;
while ~DONE_CHOOSING
  fprintf('Do you want to use discrete sample(s) to calculate offset in SUNA data?\n')
  fprintf('  <0>  No\n')
  fprintf('  <1>  Yes - calculate offset using SINGLE discrete sample concentration(s)\n')
  fprintf('  <2>  Yes - calculate offset using INTERPOLATED discrete sample concentrations\n')
  fprintf('  <3>  Yes - calculate offset using AVERAGE discrete sample concentrations\n')
  fprintf('  <9>  STOP\n')
  chc = input('  Enter choice: ');
  if isempty(chc); chc = 0; end % Default to NO
  switch chc
    case 9
      fprintf('in keyboard mode... enter "dbcont" to continue\n')
      keyboard
      DONE_CHOOSING = 0;
    case 0
      fprintf('Not calculating offset\n')
      DONE_CHOOSING = 1;
      offset = [];
    case 1 %% SINGLE VALUE
      fprintf('\n')
      fprintf('Select discrete sample to calculate offset\n')
      num = 0;
      num_idx = [];
      fprintf('       Cruise   Station     Cast  Pressure  Distance(km)  Distance(hr)\n')
      for ncast = 1:numel(casts)
        cast_name = strrep(casts{ncast},'.','_');
        cast_name = strrep(cast_name,'-','_');
        idx_bot = mtch.(cast_name).idx_bot;
        for nsamp = 1:numel(idx_bot)
          num = num + 1;
          fprintf('  <%d>  %s \t\t %s \t %d  %.1f \t %.1fkm away \t %.1fhr away\n',num,cfg.calcasts.CRUISE{idx_bot(nsamp)},cfg.calcasts.StationID{idx_bot(nsamp)},...
            cfg.calcasts.CastID(idx_bot(nsamp)), cfg.calcasts.Pressure(idx_bot(nsamp)), ...
            cfg.calcasts.distkm(idx_bot(nsamp)), cfg.calcasts.disthr(idx_bot(nsamp)));
          num_idx(num) = idx_bot(nsamp);
        end
      end
      
      select_another_sample = 1;  num_samples = 1;
      while select_another_sample
        chc_sample = input('  Enter sample choice: ');
        discrete_sample = cfg.calcasts(num_idx(chc_sample),:);
        no3var = char(discrete_sample.NO3_variable);
        if num_samples == 1
          offset = struct();
        end
        offset.number_of_samples = num_samples;
        n = num2str(num_samples);
        offset.(['info' n]) = discrete_sample;
        cast_name = strrep(char(offset.(['info' n]).StationID),'.','_');
        cast_name = strrep(cast_name,'-','_');
        offset.(['Bottle_date' n])      = datestr(nanmean(discrete_sample.dnum));
        offset.(['Bottle_datenum' n])   = nanmean(discrete_sample.dnum);
        offset.(['Bottle_' no3var n])   = discrete_sample.(no3var);
        offset.(['SUNA_NO3_uM' n])      = mtch.(cast_name).SUNA_NO3_uM;
        offset.(['SUNA_NO3_uM_TCSS' n]) = mtch.(cast_name).SUNA_NO3_uM_TCSS;
        
        select_another_sample = input('Do you want to calculate offset with another discrete sample? <1/0> ');
        num_samples = num_samples + 1;
      end
      DONE_CHOOSING = 1;
    case 2 %% INTERPOLATED VALUE
      select_another_sample = 1;  num_samples = 1;
      while select_another_sample
        fprintf('Select which cruise/cast to calculate interpolated offset\n')
        fprintf('       Cruise   Station     Cast  Pressure  Distance(km)  Distance(hr)\n')
        num = 0;
        for ncast = 1:numel(casts)
          cast_name = strrep(casts{ncast},'.','_');
          cast_name = strrep(cast_name,'-','_');
          idx_bot = mtch.(cast_name).idx_bot;
          num = num + 1;
          fprintf('  <%d>  %s \t\t %s \t %.1fkm away \t %.1fhr away\n',num,cfg.calcasts.CRUISE{idx_bot(1)},cfg.calcasts.StationID{idx_bot(1)},cfg.calcasts.distkm(idx_bot(1)), cfg.calcasts.disthr(idx_bot(1)));
        end
        chc_cruise = input('  Enter which cast: ');
        idx_bot = mtch.(strrep(casts{chc_cruise},'.','_')).idx_bot;
        fprintf('\n  Select sample range to use for interpolation \n')
        for nsamp = 1:numel(idx_bot)
          fprintf('    <N%d>  %s \t\t %s \t %d  %.1f \t %.1fkm away \t %.1fhr away\n',nsamp,cfg.calcasts.CRUISE{idx_bot(nsamp)},cfg.calcasts.StationID{idx_bot(nsamp)},...
            cfg.calcasts.CastID(idx_bot(nsamp)), cfg.calcasts.Pressure(idx_bot(nsamp)), ...
            cfg.calcasts.distkm(idx_bot(nsamp)), cfg.calcasts.disthr(idx_bot(nsamp)));
        end
        range_choice_done = 0;
        while ~range_choice_done
          chc_samples1 = input('  Enter first N#: ');
          chc_samples2 = input('  Enter last  N#: ');
          if ~isnumeric(chc_samples1) || ~isnumeric(chc_samples1) || isequal(chc_samples1,chc_samples2)
            range_choice_done = 0;
            fprintf('incorrect indice choices, try again\n')
            continue
          end
          try
            if chc_samples2 < chc_samples1
              discrete_samples = cfg.calcasts(idx_bot(chc_samples2:chc_samples1),:);
            else
              discrete_samples = cfg.calcasts(idx_bot(chc_samples1:chc_samples2),:);
            end
            range_choice_done = 1;
          catch
            fprintf('indice choices did not work, try again\n')
            range_choice_done = 0;
            continue
          end
        end
        n = num2str(num_samples);
        if num_samples == 1
          offset = struct();
        end
        no3var    = discrete_samples.NO3_variable{1};
        cast_name = discrete_samples.StationID{1};
        cast_name = strrep(cast_name,'.','_');
        cast_name = strrep(cast_name,'-','_');
        offset.(['info' n]) = discrete_samples;
        offset.number_of_samples = num_samples;
        offset.(['Bottle_date' n])      = datestr(nanmean(discrete_samples.dnum));
        offset.(['Bottle_datenum' n])   = nanmean(discrete_samples.dnum);
        offset.(['Bottle_' no3var n])   = discrete_samples.(no3var)';
        offset.(['Bottle_pressure' n])  = discrete_samples.Pressure';
        offset.(['Bottle_' no3var '_interpolated' n]) = round(interp1(discrete_samples.Pressure,discrete_samples.(no3var),mtch.(cast_name).SUNA_pressure),2);
        offset.(['SUNA_DateRange' n])   = [datestr(data_avg.datenum(mtch.(cast_name).idx_SUNA(1))) ' - ' datestr(data_avg.datenum(mtch.(cast_name).idx_SUNA(end)))];
        offset.(['SUNA_NO3_uM' n])      = mtch.(cast_name).SUNA_NO3_uM;
        offset.(['SUNA_NO3_uM_TCSS' n]) = mtch.(cast_name).SUNA_NO3_uM_TCSS;
        offset.(['SUNA_pressure' n])    = mtch.(cast_name).SUNA_pressure;
        % Loop through again or end script
        select_another_sample = input('Do you want to calculate offset with another discrete sample? <1/0> ');
        num_samples = num_samples + 1;
      end
      DONE_CHOOSING = 1;
    case 3 %% AVERAGE VALUE
      select_another_sample = 1;  num_samples = 1;
      while select_another_sample
        fprintf('Select which cruise/cast to calculate averaged offset\n')
        fprintf('       Cruise   Station     Cast  Pressure  Distance(km)  Distance(hr)\n')
        num = 0;
        for ncast = 1:numel(casts)
          cast_name = strrep(casts{ncast},'.','_');
          cast_name = strrep(cast_name,'-','_');
          idx_bot = mtch.(cast_name).idx_bot;
          num = num + 1;
          fprintf('  <%d>  %s \t\t %s \t %.1fkm away \t %.1fhr away\n',num,cfg.calcasts.CRUISE{idx_bot(1)},cfg.calcasts.StationID{idx_bot(1)},cfg.calcasts.distkm(idx_bot(1)), cfg.calcasts.disthr(idx_bot(1)));
        end
        chc_cruise = input('  Enter which cast: ');
        cast_name = strrep(casts{chc_cruise},'.','_');
        cast_name = strrep(cast_name,'-','_');
        idx_bot = mtch.(cast_name).idx_bot;
        fprintf('\n  Select sample range to use for average \n')
        for nsamp = 1:numel(idx_bot)
          fprintf('    <N%d>  %s \t\t %s \t %d  %.1f \t %.1fkm away \t %.1fhr away\n',nsamp,cfg.calcasts.CRUISE{idx_bot(nsamp)},cfg.calcasts.StationID{idx_bot(nsamp)},...
            cfg.calcasts.CastID(idx_bot(nsamp)), cfg.calcasts.Pressure(idx_bot(nsamp)), ...
            cfg.calcasts.distkm(idx_bot(nsamp)), cfg.calcasts.disthr(idx_bot(nsamp)));
        end
        range_choice_done = 0;
        while ~range_choice_done
          chc_samples1 = input('  Enter first N#: ');
          chc_samples2 = input('  Enter last  N#: ');
          if ~isnumeric(chc_samples1) || ~isnumeric(chc_samples1) || isequal(chc_samples1,chc_samples2)
            range_choice_done = 0;
            fprintf('incorrect indice choices, try again\n')
            continue
          end
          try
            if chc_samples2 < chc_samples1
              discrete_samples = cfg.calcasts(idx_bot(chc_samples2:chc_samples1),:);
            else
              discrete_samples = cfg.calcasts(idx_bot(chc_samples1:chc_samples2),:);
            end
            range_choice_done = 1;
          catch
            fprintf('indice choices did not work, try again\n')
            range_choice_done = 0;
            continue
          end
        end
        n = num2str(num_samples);
        if num_samples == 1
          offset = struct();
        end
        no3var    = discrete_samples.NO3_variable{1};
        cast_name = discrete_samples.StationID{1};
        cast_name = strrep(cast_name,'.','_');
        cast_name = strrep(cast_name,'-','_');
        offset.(['info' n]) = discrete_samples;
        offset.number_of_samples = num_samples;
        offset.(['Bottle_date' n])      = datestr(nanmean(discrete_samples.dnum));
        offset.(['Bottle_datenum' n])   = nanmean(discrete_samples.dnum);
        offset.(['Bottle_' no3var n])   = discrete_samples.(no3var)';
        offset.(['Bottle_pressure' n])  = discrete_samples.Pressure';
        offset.(['Bottle_' no3var '_averaged' n]) = round(nanmean(discrete_samples.(no3var)),2);
        offset.(['SUNA_DateRange' n])   = [datestr(data_avg.datenum(mtch.(cast_name).idx_SUNA(1))) ' - ' datestr(data_avg.datenum(mtch.(cast_name).idx_SUNA(end)))];
        offset.(['SUNA_NO3_uM' n])      = mtch.(cast_name).SUNA_NO3_uM;
        offset.(['SUNA_NO3_uM_TCSS' n]) = mtch.(cast_name).SUNA_NO3_uM_TCSS;
        offset.(['SUNA_pressure' n])    = mtch.(cast_name).SUNA_pressure;
        % Loop through again or end script
        select_another_sample = input('Do you want to calculate offset with another discrete sample? <1/0> ');
        num_samples = num_samples + 1;
      end
      DONE_CHOOSING = 1;
  end %% SWITCH CHOICE
end %% WHILE DONE_CHOOSING

% IF DECIDED NOT TO CALCULATE OFFSET, pass through empty variable
if ~isvarname('offset')
  offset = [];
end

end %% MAIN FUNCTION

