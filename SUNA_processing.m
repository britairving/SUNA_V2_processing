function [cfg,data_proc,meta_proc,data_avg,meta_avg] = SUNA_processing(cfg,data_pre,meta_pre)
%FUNCTION SUNA_processing 
%
%  Syntax:
%    [cfg,data_proc,meta_proc] = SUNA_processing(cfg,data_pre,meta_pre)
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
%% 0 | Script set up
dbstop if error

%% 0 | Limit for testing/building processing workflow
cfg.testing = 0;
if cfg.testing
  data_pre = struct2table(data_pre);
  data_pre = data_pre(1:20000,:);% data_pre(1:10:end,:);% [data_pre(1:5000,:);data_pre(end-5000:end,:)];
  data_pre = table2struct(data_pre,'ToScalar',1);
end

%% 1 | Initialize data_proc and meta_proc structures
% meta_proc will have all the same fieldnames as data_proc and store
% information about variable such as units, conversions, etc.
data_proc = data_pre; 
meta_proc = meta_pre;

%% 2 | Define flag variable
% Initialize all flag values to "not evaluated"
data_proc.flag = repmat(cfg.flag.not_evaluated,size(data_proc.datetime));
% update metadata
meta_proc.flag = struct();
meta_proc.flag = cfg.flag;

%% 3 | Discard data before deployment and after recovery
if strcmp(cfg.SUNA.ctd.level,'Level 2')
  % Theoretically, CTD data should not exist  before/after deployment so
  % just flag based on finite ctd data
  makefig; ax = gca;
  idx_fin1 = find(isfinite(data_proc.salinity),1);
  idx_fin2 = find(isfinite(data_proc.salinity),1,'last');
  yyaxis(ax,'left'); ax.YAxis(1).Color = 'k'; 
  ax.YLabel.String = 'NO_3 [uM]';
  plot(ax,data_proc.datenum,data_proc.NO3_uM,'ks','LineWidth',2);
  hold(ax,'on'); grid(ax,'on'); ax.YDir = 'normal';
  h1 = plot(ax,data_proc.datenum(1:idx_fin1),  data_proc.NO3_uM(1:idx_fin1),  'sr','MarkerFaceColor','r');
  h2 = plot(ax,data_proc.datenum(idx_fin2:end),data_proc.NO3_uM(idx_fin2:end),'sr','MarkerFaceColor','r');
  ax.YLim = [-15 40];
  text(ax,0.01,0.98,'Salinity','FontSize',16,'Color','b','Units','Normalized','FontWeight','bold','BackgroundColor','w')
  text(ax,0.01,0.94,'NO_3 [uM]', 'FontSize',16,'Color','k','Units','Normalized','FontWeight','bold','BackgroundColor','w')
  text(ax,0.01,0.90,'Before/After deployment','FontSize',16,'Color','r','Units','Normalized','FontWeight','bold','BackgroundColor','w')
  datetick(ax,'x','mmm','keeplimits');
  yyaxis(ax,'right'); ax.YAxis(2).Color  = 'b';
  plot(ax,data_proc.datenum,data_proc.salinity,'b*','LineWidth',2);
  ax.YLabel.String = 'Salinity';
  yyaxis(ax,'left');
  done = 0;
  while ~done
    fprintf('------------------------------------------------------------\n')
    fprintf('Automatically flagged data from before and after deployment are highlighted in red\n')
    fprintf('------------------------------------------------------------\n')
    fprintf('  <1> Flag all highlighted SUNA data (default)\n')
    fprintf('  <2> Manually highlight data to flag\n')
    fprintf('  <3> Do nothing\n')
    fprintf('  <9> Stop\n')
    chc_rm = input('  Enter choice: ');
    if isempty(chc_rm) % default to yes
      chc_rm = 1;
    end
    switch chc_rm
      case 1  % Flag all automatically highlighted data
        done = 1;
        data_proc.flag(1:idx_fin1)   = meta_proc.flag.bad;
        data_proc.flag(idx_fin2:end) = meta_proc.flag.bad;
      case 2 % Man
        % remove highlighted data
        delete(h1); delete(h2);
        % store original limits
        xlims = ax.XLim;
        ylims = ax.YLim;
        dnum  = data_proc.datenum;
        % zoom into first 2 days
        ax.XLim = [dnum(1)-1 dnum(1)+2];
        ax.YLim = ylims;
        done_pre = 0;
        while ~done_pre
          fprintf('-------- select range of data to remove before deployment --------\n')
          [rm_pre,~] = ginput(1);
          ix_pre = dnum <= rm_pre;
          h1 = plot(ax,data_proc.datenum(ix_pre),data_proc.NO3_uM(ix_pre),'sr','MarkerFaceColor','r');
          done_pre = input(' is that correct? <1/0> ');
        end
        % zoom into last 2 days
        ax.XLim = [dnum(end)-2 dnum(end)+1];
        ax.YLim = ylims;
        done_post = 0;
        while ~done_post
          fprintf('-------- select range of data to remove after recovery --------\n')
          [rm_post,~] = ginput(1);
          ix_post = dnum >= rm_post;
          h2 = plot(ax,data_proc.datenum(ix_post),data_proc.NO3_uM(ix_post),'sr','MarkerFaceColor','r');
          done_post = input(' is that correct? <1/0> ');
        end
        % Reset limits
        ax.XLim = xlims;
        ax.YLim = ylims;
        % flag data
        data_proc.flag(ix_pre)  = meta_proc.flag.bad;
        data_proc.flag(ix_post) = meta_proc.flag.bad;
        done = 1;
      case 3
        done = 1;
      case 9
        fprintf('in keyboard mode...\n')
        keyboard
        done = 0;
    end
  end
else
  % Set this based on salinity values less than 30 PSU
  if isfield(cfg.SUNA,'depth')
    flag_depth = cfg.SUNA.depth - round(cfg.SUNA.depth/3)*2;
  elseif isfield(data_proc,'pressure')
    p = nanmean(data_proc.pressure);
    flag_depth = p - round(p/3)*2;
  end
  
  notyet_deployed = find(data_proc.salinity < 20 | data_proc.pressure < flag_depth);
  % Set CTD and most SUNA fields to NaN
  % data_proc.NO3_uM(notyet_deployed)   = NaN;
  % data_proc.NO3_mgNL(notyet_deployed) = NaN;
  data_proc.flag(notyet_deployed) = meta_proc.flag.bad;
  % Catch others by outliers in the 1 week previous to recovery as flagged by
  % salinity. Idea is to catch when mooring is starting to be dragged up.
  after_recovery = notyet_deployed(find(diff(notyet_deployed > 1000)) + 1) - 1;
  if ~isempty(after_recovery)
    oneweek_before = find(data_proc.datenum > data_proc.datenum(after_recovery)-7 & data_proc.datenum < data_proc.datenum(after_recovery));
    after_recovery = oneweek_before(isoutlier(data_proc.pressure(oneweek_before)));
    if ~isempty(after_recovery)
      % Just go to the end of the deployment in case outliers are not flagged in
      % pressure because pockets of data are all the same
      after_recovery = after_recovery(1):numel(data_proc.pressure);
      % % check on plot
      makefig; ax = gca ;hold(ax,'on'); grid(ax,'on'); ax.YDir = 'normal';
      plot(data_proc.datenum(data_proc.flag <= meta_proc.flag.not_evaluated),data_proc.pressure(data_proc.flag <= meta_proc.flag.not_evaluated),'*');
      plot(data_proc.datenum(after_recovery),data_proc.pressure(after_recovery),'kd','MarkerFaceColor','y','LineWidth',1)
      % Set flag to bad
      data_proc.flag(after_recovery) = meta_proc.flag.bad;
    end
  end
end

%% 4 | Calculate measurement data_avg (e.g. bursts)
% Instrument (usually) operating mode is periodic. Each measurement
% sequence starts with shutter closed where "dark" records zero and then
% records counts due to thermal noise. Use this to separate measurement
% segments. 
idx_periodic = data_proc.dark == 0;
seg_start = find(idx_periodic == 1);
data_proc.burst_index = zeros(size(data_proc.datenum));
% update metadata
meta_proc.burst_index = struct();
meta_proc.burst_index.meaning = 'measurement bursts, or segments, as sequential integers';
meta_proc.burst_index.sources = 'dark';
meta_proc.burst_index.units   = 'integer';
for nseg = 1:numel(seg_start)
  % Store the shutter closed indice as segment integer + 0.5
  data_proc.burst_index(seg_start(nseg)) = nseg - 0.5;
  if nseg == numel(seg_start)
    seg_rng = seg_start(nseg)+1:numel(data_proc.datenum);
  else
    seg_rng = seg_start(nseg)+1:seg_start(nseg+1)-1;
  end
  data_proc.burst_index(seg_rng) = nseg;
end
% Update metadata to include number of segments
meta_proc.burst_index.number_of_bursts = max(data_proc.burst_index);
% Check that segments identified as expected
segs = unique(data_proc.burst_index(data_proc.dark ~= 0));
if any(diff(segs) > 1)
  % this can happen because extra files were read and include data pre/post
  % deployment
  fprintf('Segments not increasing by 1 as expected, check into this\n');
  keyboard
end

%% 5 | Automated QC
data_proc = SUNA_auto_qc(cfg,data_proc,meta_proc,cfg.SUNA.cal1.caldata);

%% 6 | Correction | Temperature compensated, salinity subtracted (TCSS)
%   Sakamoto, C. M., K. S. Johnson, and L. J. Coletti (2009), Improved
%   algorithm for the computation of nitrate concentrations in seawater
%   using an in situ ultraviolet spectrophotometer, Limnol. Oceanogr.
%   Methods, 7, 132–143.
if ~cfg.testing
  save_mid = fullfile(cfg.datadir,[cfg.project '_Sakamotoetal2009_full.mat']);
else
  save_mid = fullfile(cfg.datadir,[cfg.project '_Sakamotoetal2009_part.mat']);
end
% if exist(save_mid,'file')
%   fprintf('Loading %s\n',save_mid)
%   load(save_mid,'NTR','rmse');
% else
  tic
  fprintf('Beginning Sakamoto et al 2009 ISUS TS corrections\n')
  if cfg.SUNA.TCSS_ignore_flag  % do not calculate TCSS where data has been flagged
    [NTR,rmse] = ISUS_REPROCESSOR_v2_bi(data_proc,meta_proc,cfg.SUNA.cal1.caldata); % pass meta_proc if want to ignore bad data
  else % calculate TCSS for every measurement regardless of flag
    [NTR,rmse] = ISUS_REPROCESSOR_v2_bi(data_proc,[],cfg.SUNA.cal1.caldata);        % do NOT pass meta_proc through if want to include bad data
  end
  toc
  fprintf('saving data corrected following Sakamoto et al 2009 to file: %s\n',save_mid);
  save(save_mid,'NTR','rmse','-v7.3');
% end

% Store corrected values to data_proc structure
data_proc.NO3_uM_TCSS      = NTR;
data_proc.NO3_uM_TCSS_rmse = rmse;
% update metadata
meta_proc.NO3_uM_TCSS         = meta_proc.NO3_uM;
meta_proc.NO3_uM_TCSS.sources = {'NO3_uM' 'dark' 'salinity' 'temperature' 'spectrum_channels'};
meta_proc.NO3_uM_TCSS.method  = 'ISUS_REPROCESSOR_v2_bi.m';
meta_proc.NO3_uM_TCSS.calfile = cfg.SUNA.cal1.calfile;
meta_proc.NO3_uM_TCSS.caldata = cfg.SUNA.cal1.caldata;
meta_proc.NO3_uM_TCSS.meaning = 'Sakamoto et al. (2009) Temperature compensated salinity subtracted (TCSS)';
meta_proc.NO3_uM_TCSS_rmse          = struct();
meta_proc.NO3_uM_TCSS_rmse.sources = {'NO3_uM' 'dark' 'salinity' 'temperature' 'spectrum_channels'};
meta_proc.NO3_uM_TCSS_rmse.method  = 'ISUS_REPROCESSOR_v2_bi.m';
meta_proc.NO3_uM_TCSS_rmse.calfile = cfg.SUNA.cal1.calfile;
meta_proc.NO3_uM_TCSS_rmse.meaning = 'Sakamoto et al. (2009) TCSS RMS error of fit';

%% 7 | Average segment data
data_avg = struct();
meta_avg = struct();
data_avg.number_of_bursts = meta_proc.burst_index.number_of_bursts;
meta_avg.number_of_bursts = meta_proc.burst_index.number_of_bursts;

% initialize data_avg flag
data_avg.project = repmat(data_proc.project(1),data_avg.number_of_bursts,1);
data_avg.flag    = ones(data_avg.number_of_bursts,1)*meta_proc.flag.not_evaluated; % store flags
% initialize variables
vars = {'datenum' 'datetime' 'burst_index' 'NO3_uM' 'fit_RMSE'  'NO3_uM_TCSS' 'NO3_uM_TCSS_rmse' 'pressure' 'salinity' 'temperature' 'density' 'dark' 'cum_lampon_sec' 'cdom' 'fluor_mgm3' 'PAR' 'turbidity'};
vars = vars(ismember(vars,fieldnames(data_proc)));
% Initialize variables
for nvar = 1:numel(vars)
  vname = vars{nvar};
  if strcmp(vname,'datetime')
    data_avg.(vname) = NaT(data_avg.number_of_bursts,1);
  else
    data_avg.(vname) = nan(data_avg.number_of_bursts,1);
  end
  meta_avg.(vname)   = meta_proc.(vname);
end
% Loop through segments and pull out average counts and nitrate
for nseg = 1:data_avg.number_of_bursts
  rng = find(data_proc.burst_index == nseg & data_proc.flag <= meta_proc.flag.not_evaluated);
  for nvar = 1:numel(vars)
    vname = vars{nvar};
    if strcmp(vname,'burst_index')
      data_avg.burst_index(nseg) = nseg;
    else
      data_avg.(vname)(nseg) = mean(data_proc.(vname)(rng));
    end
  end
  if isempty(rng)
    % even if all data are flagged, calculate mean of key variables
    data_avg.datenum(nseg)        = mean(data_proc.datenum(data_proc.burst_index == nseg));
    data_avg.datetime(nseg)       = mean(data_proc.datetime(data_proc.burst_index == nseg));
    data_avg.cum_lampon_sec(nseg) = round(mean(data_proc.cum_lampon_sec(data_proc.burst_index == nseg)));
    data_avg.dark(nseg)           = round(mean(data_proc.cum_lampon_sec(data_proc.dark == nseg)));
    data_avg.flag(nseg)           = meta_proc.flag.bad;
  end
end

%% 8 | Filter and smooth data
hr_window = 35;    % 35-hour filter window for movmedian filter
Tc = hr_window*60; % minutes in 35 hours, used for Lanczos filter
% Find segment sampling interval in minutes
dT = nan(1,data_avg.number_of_bursts);
for nseg = 2:data_avg.number_of_bursts-1
  dT(nseg) = round(etime(datevec(data_avg.datenum(nseg)),datevec(data_avg.datenum(nseg-1)))./60); % Seconds to minutes
end
dT = nanmedian(dT);

% Instead.. try movmdedian with 35 hour filter
filter_vars = {'NO3_uM' 'NO3_uM_TCSS'};
idx = isfinite(data_avg.datetime);
for nv = 1:numel(filter_vars)
  smo_var  = [filter_vars{nv} '_smo'];
  % Initialize
  data_avg.(smo_var)      = nan(size(data_avg.(filter_vars{nv})));
  data_avg.(smo_var)(idx) = smoothdata(data_avg.(filter_vars{nv})(idx),'movmedian',hours(hr_window),'SamplePoints',data_avg.datetime(idx));
  data_proc.(smo_var)     = smoothdata(data_proc.(filter_vars{nv}),'movmedian',hours(hr_window),'SamplePoints',data_proc.datetime);
  % Update metadata
  meta_proc.(smo_var)  = meta_proc.(filter_vars{nv});
  meta_proc.(smo_var).sources = filter_vars{nv};
  meta_proc.(smo_var).method  = 'smoothdata.m, movmedian';
  meta_proc.(smo_var).window  = [num2str(hr_window) 'hours'];
  % Update metadata
  meta_avg.(smo_var)  = meta_avg.(filter_vars{nv});
  meta_avg.(smo_var).sources = filter_vars{nv};
  meta_avg.(smo_var).method  = 'smoothdata.m, movmedian';
  meta_avg.(smo_var).window  = [num2str(hr_window) 'hours'];
  
%   % Lanczos Filter [Carlos Adrian Vargas Aguilera (2020). LanczosFilter.m
%   % (https://www.mathworks.com/matlabcentral/fileexchange/14041-lanczosfilter-m), MATLAB Central File Exchange. Retrieved October 26, 2020. ]
%   % Nearest could find to Mordy et al's [2020] 35 h, cosine-squared,
%   % tapered Lanczos filter to remove tidal and higher-frequency
%   % variability.
%   smo_var2 = [filter_vars{nv} '_smo2'];
%   data_avg.(smo_var2)      = nan(size(data_avg.(filter_vars{nv})));
%   data_avg.(smo_var2)(idx) = lanczosfilter(data_avg.(filter_vars{nv})(idx),dT,1/Tc,[],'low');
%   data_proc.(smo_var2)     = lanczosfilter(data_proc.(filter_vars{nv}),dT,1/Tc,[],'low');
%   % Update metadata
%   meta_proc.(smo_var)  = meta_proc.(filter_vars{nv});
%   meta_proc.(smo_var).sources = filter_vars{nv};
%   meta_proc.(smo_var).method  = 'smoothdata.m, movmedian';
%   meta_proc.(smo_var).window  = [num2str(hr_window) 'hours'];
%   
%   meta_proc.(smo_var2) = meta_proc.(filter_vars{nv});
%   meta_proc.(smo_var2).sources = filter_vars{nv};
%   meta_proc.(smo_var2).method  = 'lanczosfilter.m';
%   meta_proc.(smo_var2).reference = 'Carlos Adrian Vargas Aguilera (2020). LanczosFilter.m (https://www.mathworks.com/matlabcentral/fileexchange/14041-lanczosfilter-m), MATLAB Central File Exchange. Retrieved October 26, 2020.';
%   meta_proc.(smo_var2).dT   = [num2str(dT) 'minutes, Sampling interval'];
%   meta_proc.(smo_var2).Cf   = ['1/Tc, Tc=' num2str(hr_window) '(hours)*60(minutes/hours), Cut-off frequency'];
%   meta_proc.(smo_var2).M    = '[], Number of coefficients (default: 100)';
%   meta_proc.(smo_var2).pass = 'low-pass';
end

end %% function [cfg, t] = SeapHOx_find_zero(cfg,data)