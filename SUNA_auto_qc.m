function [data] = SUNA_auto_qc(cfg,data_proc,meta_proc,cal)
%FUNCTION SUNA_processing
%
%  Syntax:
%    [cfg,data_proc,meta_proc] = SUNA_processing(cfg,data,meta_pre,data_raw,meta_raw)
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
%   %INTERFERENCE
%   SUNA V2 MANUAL pg 18
%   An interfering species generates a spurious nitrate concentration when the spectral
%   characteristics of the interfering species resembles that of nitrate. Typically, an RMSE
%   value that is more than a few times the RMSE of a pure nitrate sample should be taken
%   as an indication that interfering species are impacting the measurement. The RMSE
%   value is the square root of the mean of the sum of the squared differences between the
%   measured and the fitted absorbance; it provides a measure for the quality of the fit.
%   Independent measurements of turbidity and CDOM, as well as an analysis of the
%   absorption spectrum, can refine the impact analysis.
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%

%% 0 | Script set up
dbstop if error
close all
auto_qc = 1; % 0 = asks user to participate, 1 = runs all throughs through without user interrupts
if isfield(cfg,'save_figures')
  savefig = cfg.save_figures;
else
  savefig = false; % default to no
end

segment_clr = 'y';
if strcmp(cfg.project,'CEO_2016')
  percent_threshold = -7; % percent difference threshold [%]
elseif strcmp(cfg.project,'CEO_2018')
  percent_threshold = -7; % percent difference threshold [%]
else
  percent_threshold = -10; % percent difference threshold [%]
end
std_threshold = 1.5; % Used for test #5 (Outliers based on deviations of percent_difference_threshold between points and meeting the criteria > avg + std_threshold*std\)
%% Initialize data structure so do not overwrite any data_proc variables
data = data_proc;
% definte indicies for key wavelengths
i254_and_350  = find(round(cal.Wavelength)==254 | round(cal.Wavelength)==350);
i217_to_240   = find(cal.Wavelength >= 217 & cal.Wavelength < 241);

%% Calculate average spectra for each measurement segment/burst 
calculate_segment_averages

%% QC Test #1 | Absorbance at 254nm or 350nm > 1.3
a = plot_fullresolution_3axes(i254_and_350); % plot density, NO3, spectrum channel, percent differences and return axes
a(3).Title.String = {'QC test: #1 | Absorbance at 254nm or 350nm > 1.3'};
% Absorption: The data output of the SUNA V2 is the absorption at 350 nm
% and 254 nm (A350 and A254). These wavelengths are outside the nitrate
% absorption range and can be used to make an estimate of the impact of
% CDOM. If absorption is high (>1.3 AU), the SUNA will not be able to
% collect adequate light to make a measurement. [ SUNA V2 MANUAL ]
bad_abs = data.abs_254nm > 1.3 | data.abs_350nm > 1.3;
plot(a(2),data.datenum(bad_abs),data.NO3_uM(bad_abs),'g*','MarkerSize',4,'linewidth',2);
plot(a(3),data.datenum(bad_abs),data.spectrum_channels(bad_abs,i254_and_350),'g*','MarkerSize',4,'linewidth',2);
data.flag(bad_abs) = meta_proc.flag.bad;
% flag segments as questionable
questionable_segments = unique(data.burst_index(bad_abs));
bad_abs_segments = ismember(data.burst_index,questionable_segments) & data.flag <= meta_proc.flag.not_evaluated; % Don't overwrite/plot the already flagged values
data.flag(bad_abs_segments) = meta_proc.flag.questionable;
data_avg.flag(questionable_segments) = meta_proc.flag.questionable;
if savefig
  savename = fullfile(cfg.datadir,[cfg.project '_qctest_abs254_abs350']);
  standard_printfig_highrespng(savename);
end

%% QC Test #2 | Nitrate spectral fit RMSE > 1e-3
a = plot_fullresolution_3axes(i254_and_350); % plot density, NO3, spectrum channel, percent differences and return axes
a(3).Title.String = {'QC test: #2 |  Nitrate spectral fit RMSE > 1e-3'};
% RMSE: The root mean square error parameter from the SUNA V2 can be used
% to make an estimate of how well the nitrate spectral fit is. This should
% usually be less than 1E-3. If it is higher, there is spectral shape
% (likely due to CDOM) that adversely impacts the nitrate estimate. 
% [ SUNA V2 MANUAL ]
bad_fit = data.fit_RMSE > 1e-3 & data.flag <= meta_proc.flag.not_evaluated; % Don't overwrite/plot the already flagged values
plot(a(2),data.datenum(bad_fit),data.NO3_uM(bad_fit),'g*','MarkerSize',4,'linewidth',2);
plot(a(3),data.datenum(bad_fit),data.spectrum_channels(bad_fit,i254_and_350),'g*','MarkerSize',4,'linewidth',2);
data.flag(bad_fit) = meta_proc.flag.bad;
% flag segments as questionable
questionable_segments = unique(data.burst_index(bad_fit));
bad_fit_segments = ismember(data.burst_index,questionable_segments) & data.flag <= meta_proc.flag.not_evaluated; % Don't overwrite/plot the already flagged values
data.flag(bad_fit_segments) = meta_proc.flag.questionable;
data_avg.flag(questionable_segments) = meta_proc.flag.questionable;
if savefig
  savename = fullfile(cfg.datadir,[cfg.project '_qctest_fitRMSE']);
  standard_printfig_highrespng(savename);
end

%% QC Test #3 | Segment dark count spikes
fprintf('--------------------------------------------------------\n')
fprintf('QC test: #3 | Segment dark count spikes\n')
fprintf('--------------------------------------------------------\n')
bad_dark_segments = [];
for nseg = 1:data_avg.number_of_bursts
  % Flag all segments where dark counts were not the same
  rng_seg = find(data.burst_index == nseg);
  if ~all(data.dark(rng_seg) == data.dark(rng_seg(1)))
    data_avg.flag(nseg) = meta_proc.flag.bad;
    bad_dark_segments = [bad_dark_segments; nseg];
  end
end
if ~isempty(bad_dark_segments)
  idx_bad_dark_segments = ismember(data.burst_index,bad_dark_segments) & data.flag <= meta_proc.flag.not_evaluated;
  a = plot_fullresolution_3axes(i254_and_350); % plot density, NO3, spectrum channel, percent differences and return axes
  a(3).Title.String = {'QC test: #3 | Segment dark count spikes'};
  plot(a(2),data.datenum(idx_bad_dark_segments),data.NO3_uM(idx_bad_dark_segments),'g*','MarkerSize',4,'linewidth',2);
  plot(a(3),data.datenum(idx_bad_dark_segments),data.spectrum_channels(idx_bad_dark_segments,i254_and_350),'g*','MarkerSize',4,'linewidth',2);
  % set to bad flag
  data.flag(idx_bad_dark_segments) = meta_proc.flag.bad;
  if savefig
    savename = fullfile(cfg.datadir,[cfg.project '_qctest_darkspikes']);
    standard_printfig_highrespng(savename);
  end
end

%% QC Test #4 | Deviations in segments at 254nm or 350nm based on percent differences
fprintf('--------------------------------------------------------\n')
fprintf('QC test: #4 | Deviations of %d%% between segments\n',percent_threshold)
fprintf('--------------------------------------------------------\n')
% recalculate segments so flagged data are excluded
calculate_segment_averages
% Now calculate percent differences
calculate_percent_differences(percent_threshold); % Default percent_threshold = -10%
% Plot
a = plot_segment_percent_diff(i254_and_350); % plots (full resolution & segment) density, NO3, spectrum channel, percent differences and return axes
str = a(4).Title.String;
a(4).Title.String = {['QC test: #4 | Deviations of ' num2str(pdiff.percent_diff_threshold) '% between segments'];str};

% Loop through segments and find bad segments
number_bad_segments   = []; % Store number of flagged segments so can stop endless loop that seems to catch okay data segmemts
find_segment_outliers = 1;
while find_segment_outliers
  % Identify bad segments based on 254nm or 350nm percent differences
  % relative to other segments
  idx_bad_seg = find(any(pdiff.spec_seg(:,i254_and_350) <= pdiff.percent_diff_threshold,2));
  % find bad segments in full resoltuion
  idx_bad_segments = ismember(data.burst_index, idx_bad_seg);

  if numel(idx_bad_seg) == 0
    fprintf('No bad segments identified, exiting loop\n');
    find_segment_outliers = 0;
    continue
  else
    number_bad_segments = [number_bad_segments; numel(idx_bad_seg)];
  end
  % Catch if looping through endlessly and flagging good data
  % happened with CEO 2017 
  if numel(number_bad_segments) > 3 && all(numel(idx_bad_seg) == number_bad_segments(end-3:end))
    fprintf('Maybe caught in repeated cycle of flagging okay data... exiting loop\n');
    find_segment_outliers = 0;
    continue
  end
  % show bad segments
  plot(a(1),data.datenum(idx_bad_segments),data.density(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
  plot(a(2),data.datenum(idx_bad_segments),data.NO3_uM(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
  plot(a(3),data.datenum(idx_bad_segments),data.spectrum_channels(idx_bad_segments,i254_and_350),'g*','MarkerSize',3,'linewidth',2);
  plot(a(4),data.datenum(idx_bad_segments),pdiff.spec_full(idx_bad_segments,i254_and_350),'g*','MarkerSize',3,'linewidth',2);
  
  % Flag bad segments in full resoltution, and segment data
  data_avg.flag(idx_bad_seg)  = meta_proc.flag.bad;
  data.flag(idx_bad_segments) = meta_proc.flag.bad;
 
  % recalculate percent differences
  calculate_percent_differences(pdiff.percent_diff_threshold); % Default percent_threshold = -10%
  
  % Replot percent differences on top axes
  cla(a(4))
  plot(a(4),data.datenum(pdiff.idx_ful_good),pdiff.spec_full(pdiff.idx_ful_good,i254_and_350),'*','MarkerSize',3,'LineWidth',2);
  hold(a(4),'on'); grid(a(4),'on'); ylim(a(4),[-100 100])
  plot(a(4),data_avg.datenum,pdiff.spec_seg(:,i254_and_350),'*','Color',segment_clr,'MarkerSize',2,'lineWidth',1);
  fprintf('\n')
  fprintf('Percent difference threshold: %d\n',pdiff.percent_diff_threshold);
  fprintf('Number of bad segments: %d\n',numel(idx_bad_seg));
  chc_done = 0;
  if ~auto_qc
    while ~chc_done
      fprintf('  <1> Loop through and find more bad segments\n')
      fprintf('  <0> Exit loop\n')
      fprintf('  <9> Enter keyboard mode\n')
      chc_done = input('  Enter choice: ');
      if isempty(chc_done); chc_done = 1; end % default to find more outliers
      switch chc_done
        case 0
          find_segment_outliers = 0;
          chc_done = 1;
        case 1
          find_segment_outliers = 1;
          chc_done = 1;
        case 9
          keyboard
          chc_done = 0;
      end
    end
  end
end
if savefig
  savename = fullfile(cfg.datadir,[cfg.project '_qctest_segment_outliers_at_254nm_350nm']);
  standard_printfig_highrespng(savename);
end

%% QC Test #5 | Single Point deviations at 254nm or 350nm based on percent differences
fprintf('--------------------------------------------------------\n')
fprintf('QC test: #5 | Outliers based on deviations of %d%% between points and meeting the criteria > avg + %.1f*std\n',percent_threshold,std_threshold)
fprintf('--------------------------------------------------------\n')
% "A data point was discarded if its standard deviation (based on 13
% samples for each burst) exceeded 0.2 ?M, which was coincident with low
% detector intensities." REF: [Randelhoff, A., Sundfjord, A. & Reigstad, M.
% Seasonal variability and fluxes of nitrate in the surface waters over the
% Arctic shelf slope. Geophys. Res. Lett. 42, 3442–3449 (2015).]
% Loop through segments and find bad segments
a = plot_fullresolution_percent_diff(i254_and_350); % plot density, NO3, spectrum channel, percent differences and return axes
str = a(4).Title.String;
a(4).Title.String = {['QC test: #5 | Outliers based on deviations of ' num2str(percent_threshold) '% between points and meeting the criteria > avg + ' num2str(std_threshold) '*std'];str};
find_segment_outliers = 1;
while find_segment_outliers
   idx_bad_pts = find(any(pdiff.spec_full(:,i254_and_350) <= percent_threshold,2));
  % Now make make sure points are with
  seg_idx = data.burst_index(idx_bad_pts);
  idx_bad_pts_2 = data.NO3_uM(idx_bad_pts) > data_avg.NO3_uM(seg_idx)+std_threshold*data_avg.NO3_uM_std(seg_idx) | data.NO3_uM(idx_bad_pts) < data_avg.NO3_uM(seg_idx)-std_threshold*data_avg.NO3_uM_std(seg_idx);
  % Reindex and only include those points that meet the above criteria
  idx_bad_pts   = idx_bad_pts(idx_bad_pts_2);
  if numel(idx_bad_pts) == 0
    find_segment_outliers = 0;
    continue
  else
    plot(a(2),data.datenum(idx_bad_pts),data.NO3_uM(idx_bad_pts),'g*','MarkerSize',3,'linewidth',2);
    plot(a(3),data.datenum(idx_bad_pts),data.spectrum_channels(idx_bad_pts,i254_and_350),'g*','MarkerSize',3,'linewidth',2);
    plot(a(4),data.datenum(idx_bad_pts),pdiff.spec_full(idx_bad_pts,i254_and_350),'g*','MarkerSize',3,'linewidth',2);
    % Store flags
    data.flag(idx_bad_pts) = meta_proc.flag.bad;
    % recalculate percent differencesv
    calculate_percent_differences(percent_threshold); % Default percent_threshold = -10%
    % replot %difference axes
    cla(a(4))
    plot(a(4),data.datenum(pdiff.idx_ful_good),pdiff.spec_full(pdiff.idx_ful_good,i254_and_350),'*','MarkerSize',3,'LineWidth',2);
    hold(a(4),'on'); grid(a(4),'on'); ylim(a(4),[-100 100])
    plot(a(4),data_avg.datenum,pdiff.spec_seg(:,i254_and_350),'*','Color',segment_clr,'MarkerSize',2,'lineWidth',1);
  end
end
if savefig
  savename = fullfile(cfg.datadir,[cfg.project '_qctest_point_outliers_at_254nm_350nm']);
  standard_printfig_highrespng(savename);
end

% 
% %% QC Test #6 | Deviations in segments at 254nm or 350nm based on percent differences
% percent_threshold = -2000;
% fprintf('--------------------------------------------------------\n')
% fprintf('QC test: #6 | Absorbance 254 or 350nm : Deviations of %d%% between segments\n',percent_threshold)
% fprintf('--------------------------------------------------------\n')
% calculate_percent_differences(percent_threshold); % Default percent_threshold = -10%
% % Plot
% a = plot_segment_percent_diff_abs; % plots (full resolution & segment) density, NO3, spectrum channel, percent differences and return axes
% str = a(4).Title.String;
% a(4).Title.String = {['QC test: #4 | Deviations of ' num2str(pdiff.percent_diff_threshold) '% between segments'];str};
% 
% % Loop through segments and find bad segments
% find_segment_outliers = 1;
% while find_segment_outliers
%   % Identify bad segments based on 254nm or 350nm percent differences
%   % relative to other segments
%   idx_bad_seg = find(pdiff.abs_254nm_seg <= pdiff.percent_diff_threshold | pdiff.abs_350nm_seg  <= pdiff.percent_diff_threshold);
%   % find bad segments in full resoltuion
%   idx_bad_segments = ismember(data.burst_index, idx_bad_seg);
% 
%   if numel(idx_bad_seg) == 0
%     fprintf('No bad segments identified, exiting loop\n');
%     find_segment_outliers = 0;
%     continue
%   end
%   % show bad segments
%   plot(a(2),data.datenum(idx_bad_segments),data.NO3_uM(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
%   plot(a(3),data.datenum(idx_bad_segments),data.abs_254nm(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
%   plot(a(3),data.datenum(idx_bad_segments),data.abs_350nm(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
%   plot(a(4),data.datenum(idx_bad_segments),pdiff.abs_254nm_full(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
%   plot(a(4),data.datenum(idx_bad_segments),pdiff.abs_350nm_full(idx_bad_segments),'g*','MarkerSize',3,'linewidth',2);
%   % Flag bad segments in full resoltution, and segment data
%   keyboard
%   data_avg.flag(idx_bad_seg)  = meta_proc.flag.bad;
%   data.flag(idx_bad_segments) = meta_proc.flag.bad;
%  
%   % recalculate percent differences
%   calculate_percent_differences(pdiff.percent_diff_threshold); % Default percent_threshold = -10%
%   
%   % Replot percent differences on top axes
%   cla(a(4))
%   plot(a(4),data.datenum(pdiff.idx_ful_good),pdiff.abs_254nm(pdiff.idx_ful_good),'*','MarkerSize',3,'LineWidth',2);
%   hold(a(4),'on'); grid(a(4),'on'); ylim(a(4),[-100 100])
%   plot(a(4),data.datenum(pdiff.idx_ful_good),pdiff.abs_350nm(pdiff.idx_ful_good),'*','MarkerSize',3,'LineWidth',2);
%   plot(a(4),data_avg.datenum,pdiff.abs_254nm,'*','Color',segment_clr,'MarkerSize',2,'lineWidth',1);
%   plot(a(4),data_avg.datenum,pdiff.abs_350nm,'*','Color',segment_clr,'MarkerSize',2,'lineWidth',1);
% 
%   fprintf('\n')
%   fprintf('Percent difference threshold: %d\n',pdiff.percent_diff_threshold);
%   fprintf('Number of bad segments: %d\n',numel(idx_bad_seg));
%   chc_done = 0;
%   %if ~auto_qc
%     while ~chc_done
%       fprintf('  <1> Loop through and find more bad segments\n')
%       fprintf('  <0> Exit loop\n')
%       fprintf('  <9> Enter keyboard mode\n')
%       chc_done = input('  Enter choice: ');
%       switch chc_done
%         case 0
%           find_segment_outliers = 0;
%           chc_done = 1;
%         case 1
%           find_segment_outliers = 1;
%           chc_done = 1;
%         case 9
%           keyboard
%           chc_done = 0;
%       end
%     end
%   %end
% end
% if savefig
%   savename = fullfile(cfg.datadir,[cfg.project '_qctest2']);
%   standard_printfig_highrespng(savename);
% end
% 

%% function calculate_segment_averages
  function calculate_segment_averages
    if ~exist('segments','var')
      segments = struct();
      data_avg.number_of_bursts = meta_proc.burst_index.number_of_bursts;
      % initialize segments flag
      data_avg.flag = ones(data_avg.number_of_bursts,1)*meta_proc.flag.not_evaluated; % store flags
    end
    % initialize variables
    vars = {'datenum' 'spec_channels' 'NO3_uM' 'NO3_uM_std' 'abs_254nm' 'abs_350nm' 'fit_RMSE' 'dark' 'density'};
    % Initialize variables
    for nvar = 1:numel(vars)
      vname = vars{nvar};
      if strcmp(vname,'spec_channels')
        data_avg.(vname) = nan(data_avg.number_of_bursts,numel(cal.Wavelength));
      else
        data_avg.(vname) = nan(data_avg.number_of_bursts,1);
      end
    end
    % Loop through segments and pull out average counts and nitrate
    for nseg = 1:data_avg.number_of_bursts
      rng = find(data_proc.burst_index == nseg & data_proc.flag <= meta_proc.flag.not_evaluated);
      for nvar = 1:numel(vars)
        vname = vars{nvar};
        if strcmp(vname,'spec_channels')
          data_avg.(vname)(nseg,:) = mean(data.spectrum_channels(rng,:));
        elseif strcmp(vname,'NO3_uM_std')
          data_avg.(vname)(nseg)   = std(data_proc.NO3_uM(rng));
        else
          data_avg.(vname)(nseg)   = mean(data_proc.(vname)(rng));
        end
      end
    end
  end %% function segments = calculate_segment_averages

%% function calculate_percent_differences(percent_dif)
  function calculate_percent_differences(percent_dif)
    fprintf(' Calculating percent differences in spectrum\n')
    pdiff = struct();
    % Default percent difference threshold to -10%
    if nargin == 1
      pdiff.percent_diff_threshold = percent_dif;
    else
      pdiff.percent_diff_threshold = -10;
    end
    % Calculate percent difference based nondark values and unflagged data
    pdiff.idx_ful_good = find(data.dark ~= 0 & data.flag <= meta_proc.flag.not_evaluated);%find(data.dark ~= 0 & all(isfinite(data.spectrum_channels),2));
    pdiff.idx_seg_good = find(data_avg.flag <= meta_proc.flag.not_evaluated);
    
    % Initialize percent difference fields
    pdiff.spec_full = nan(size(data.spectrum_channels));  % Full resolution percent differences at specified wavelength
    pdiff.spec_seg  = nan(size(data_avg.spec_channels));  % Segment percent differences at specified wavelength
    pdiff.abs_254nm_seg  = nan(size(data_avg.abs_254nm)); % Segment percent differences for absorbance at 254nm
    pdiff.abs_350nm_seg  = nan(size(data_avg.abs_254nm)); % Segment percent differences for absorbance at 350nm
    pdiff.abs_254nm_full = nan(size(data.abs_254nm));     % Full resolution percent differences for absorbance at 254nm
    pdiff.abs_350nm_full = nan(size(data.abs_254nm));     % Full resolution percent differences for absorbance at 350nm
    % Calculate percent differences in absorbances
    a254_ful = data.abs_254nm(pdiff.idx_ful_good);
    a350_ful = data.abs_350nm(pdiff.idx_ful_good);
    a254_seg = data_avg.abs_254nm(pdiff.idx_seg_good);
    a350_seg = data_avg.abs_350nm(pdiff.idx_seg_good);
    
    tmp_seg_254 = round( ( (a254_seg(2:end) - a254_seg(1:end-1)) ./a254_seg(1:end-1) ) *100);
    tmp_seg_350 = round( ( (a350_seg(2:end) - a350_seg(1:end-1)) ./a350_seg(1:end-1) ) *100);
    tmp_ful_254 = round( ( (a254_ful(2:end) - a254_ful(1:end-1)) ./a254_ful(1:end-1) ) *100);
    tmp_ful_350 = round( ( (a350_ful(2:end) - a350_ful(1:end-1)) ./a350_ful(1:end-1) ) *100);
    
    pdiff.abs_254nm_seg(pdiff.idx_seg_good(2:end))  = tmp_seg_254;
    pdiff.abs_350nm_seg(pdiff.idx_seg_good(2:end))  = tmp_seg_350;
    pdiff.abs_254nm_full(pdiff.idx_ful_good(2:end)) = tmp_ful_254;
    pdiff.abs_350nm_full(pdiff.idx_ful_good(2:end)) = tmp_ful_350;
    % Loop through all wavelengths and calculate percent differences
    for nwv = 1:numel(cal.Wavelength)
      %swv = num2str(cal.Wavelength(nwv)); % string wavelength
      %fprintf(' Calculating percent differences in spectrum @ %snm\n',swv)
      % calculate percent difference
      ful_spec = data.spectrum_channels(pdiff.idx_ful_good,nwv);
      seg_spec = data_avg.spec_channels(pdiff.idx_seg_good,nwv);
      
      % relative difference = [x(2) - x(1)] ./ x(1);
      % percent difference  = relative difference * 100%
      tmp_ful  = round( ( (ful_spec(2:end) - ful_spec(1:end-1)) ./ful_spec(1:end-1) ) *100);
      tmp_seg  = round( ( (seg_spec(2:end) - seg_spec(1:end-1)) ./seg_spec(1:end-1) ) *100);
      % reindex
      pdiff.spec_full(pdiff.idx_ful_good(2:end),nwv) = tmp_ful;
      pdiff.spec_seg(pdiff.idx_seg_good(2:end), nwv) = tmp_seg;
    end
  end %% function calculate_percent_differences

%% function a = plot_fullresolution_3axes(idx_wavelengths)
  function a = plot_fullresolution_3axes(idx_wavelengths)
    % default to 217 to 240
    if nargin == 0
      idx_wavelengths = i217_to_240;%1:numel(cal.Wavelength);
    end
    switch numel(idx_wavelengths)
      case 1
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm'];
      case 2
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm & ' num2str(cal.Wavelength(idx_wavelengths(end)),'%.1f') 'nm'];
      otherwise
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm to ' num2str(cal.Wavelength(idx_wavelengths(end)),'%.1f') 'nm'];
    end
    
    % plot density, NO3, spectrum channel and return axes
    a = makefig_subplots(1,3);% spectra
    xlimits = [round(data.datenum(1)-1) round(data.datenum(end)+1)];
    idx_good = find(data.dark ~= 0 & data.flag <= meta_proc.flag.not_evaluated);
    % ------ DENSITY
    plot(a(1),data.datenum,data.density,'k-*','MarkerSize',2);
    hold(a(1),'on'); grid(a(1),'on');
    a(1).YLim = round(prctile(data.density(idx_good),[0.01 99.999])); %ylim(a(1),[1023 1027])
    ylabel(a(1),'Density'); datetick(a(1),'Keepticks');
    
    % ------ RAW NO3 [MICROMOLAR]
    plot(a(2),data.datenum(idx_good),data.NO3_uM(idx_good),'k.');
    hold(a(2),'on'); grid(a(2),'on');ylabel(a(2),'NO_3 [\muM]');
    a(2).YLim = round(prctile(data.NO3_uM(idx_good),[0.02 99.99])); %ylim(a(2),no3_limits)
    
    % ------ SPECTRUM CHANNELS
    plot(a(3),data.datenum(idx_good),data.spectrum_channels(idx_good,idx_wavelengths),'*','MarkerSize',2,'LineWidth',2);
    hold(a(3),'on'); grid(a(3),'on'); ylabel(a(3),'spectral counts');
    title(a(3),title_string),
    limits = round(prctile(data.spectrum_channels(idx_good,idx_wavelengths),[0.02 99.99]));
    a(3).YLim = [min(min(limits)) max(max(limits))]; 

    a(1).XLim = xlimits;
    a(2).XLim = xlimits;
    a(3).XLim = xlimits;
    a(2).XTickLabel = [];
    a(3).XTickLabel = [];
    % link axes so when zoom all axes will change limits
    linkaxes([a(1) a(2) a(3)],'x')
  end %% function a = plot_fullresolution_3axes

%% function a = plot_fullresolution_percent_diff(idx_wavelengths)
  function a = plot_fullresolution_percent_diff(idx_wavelengths)
    % default to 217 to 240
    if nargin == 0
      idx_wavelengths = i217_to_240;%1:numel(cal.Wavelength);
    end
    switch numel(idx_wavelengths)
      case 1
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm'];
      case 2
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm & ' num2str(cal.Wavelength(idx_wavelengths(end)),'%.1f') 'nm'];
      otherwise
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm to ' num2str(cal.Wavelength(idx_wavelengths(end)),'%.1f') 'nm'];
    end
    
    % plot density, NO3, spectrum channel, percent differences and return axes
    a = makefig_subplots(1,4);% spectra
    xlimits = [data.datenum(1)-1 data.datenum(end)+1];
    idx_good = find(data.dark ~= 0 & data.flag <= meta_proc.flag.not_evaluated);
    
    % ------ DENSITY
    plot(a(1),data.datenum,data.density,'k-*','MarkerSize',2); % plot(a(1),data.datenum(idx_good),data.density(idx_good),'k-*','MarkerSize',2);
    hold(a(1),'on'); grid(a(1),'on');
    a(1).YLim = round(prctile(data.density(idx_good),[0.01 99.999]));%ylim(a(1),[1023 1027])
    ylabel(a(1),'Density'); datetick(a(1),'Keepticks');
    
    % ------ RAW NO3 [MICROMOLAR]
    plot(a(2),data.datenum(idx_good),data.NO3_uM(idx_good),'k.');
    hold(a(2),'on'); grid(a(2),'on');ylabel(a(2),'NO_3 [\muM]');
    a(2).YLim = round(prctile(data.NO3_uM(idx_good),[0.02 99.99])); %ylim(a(2),no3_limits)
    
    % ------ SPECTRUM CHANNELS
    plot(a(3),data.datenum(idx_good),data.spectrum_channels(idx_good,idx_wavelengths),'*','MarkerSize',2,'LineWidth',2);
    hold(a(3),'on'); grid(a(3),'on'); ylabel(a(3),'spectral counts');
    
    % ------ PERCENT DIFFERENCE IN SPECTRUM CHANNELS
    plot(a(4),data.datenum(idx_good),pdiff.spec_full(idx_good,idx_wavelengths),'*','MarkerSize',2,'LineWidth',2);
    hold(a(4),'on'); grid(a(4),'on'); ylim(a(4),[-100 100])
    title(a(4),title_string),
    ylabel(a(4),'% Difference')
    
    a(1).XLim = xlimits;
    a(2).XLim = xlimits;
    a(3).XLim = xlimits;
    a(4).XLim = xlimits;
    a(2).XTickLabel = [];
    a(3).XTickLabel = [];
    a(4).XTickLabel = [];
    % link axes so when zoom all axes will change limits
    linkaxes([a(1) a(2) a(3) a(4)],'x')
  end %% function a = plot_fullresolution_percent_diff

%% function a = plot_segment_percent_diff(idx_wavelengths)
  function a = plot_segment_percent_diff(idx_wavelengths)
    % default to 217 to 240
    if nargin == 0
      idx_wavelengths = i217_to_240;%1:numel(cal.Wavelength);
    end
    switch numel(idx_wavelengths)
      case 1
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm'];
      case 2
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm & ' num2str(cal.Wavelength(idx_wavelengths(end)),'%.1f') 'nm'];
      otherwise
        title_string = ['Wavelengths: ' num2str(cal.Wavelength(idx_wavelengths(1)),'%.1f') 'nm to ' num2str(cal.Wavelength(idx_wavelengths(end)),'%.1f') 'nm'];
    end
    
    % plot density, NO3, spectrum channel, percent differences and return axes
    a = makefig_subplots(1,4);% spectra
    xlimits = [data.datenum(1)-1 data.datenum(end)+1];
    idx_ful_good = find(data.dark ~= 0 & data.flag <= meta_proc.flag.not_evaluated); % & all(isfinite(data.spectrum_channels),2));
    idx_seg_good = find(data_avg.flag <= meta_proc.flag.not_evaluated);              % & all(isfinite(data_avg.spec_channels),2));
    %% Plot full resolution data
    % ------ DENSITY
    plot(a(1),data.datenum,data.density,'k-*','MarkerSize',2); % plot(a(1),data.datenum(idx_ful_good),data.density(idx_ful_good),'k-*','MarkerSize',2);
    hold(a(1),'on'); grid(a(1),'on');
    a(1).YLim = round(prctile(data.density(idx_ful_good),[0.02 99.99])); %ylim(a(1),[1023 1027])
    ylabel(a(1),'Density'); datetick(a(1),'Keepticks');
    
    % ------ RAW NO3 [MICROMOLAR]
    plot(a(2),data.datenum(idx_ful_good),data.NO3_uM(idx_ful_good),'k*','MarkerSize',3,'LineWidth',2);
    hold(a(2),'on'); grid(a(2),'on');ylabel(a(2),'NO_3 [\muM]');
    a(2).YLim = round(prctile(data.NO3_uM(idx_ful_good),[0.02 99.99])); %ylim(a(2),no3_limits) 
    
    % ------ SPECTRUM CHANNELS
    plot(a(3),data.datenum(idx_ful_good),data.spectrum_channels(idx_ful_good,idx_wavelengths),'*','MarkerSize',3,'LineWidth',2);
    hold(a(3),'on'); grid(a(3),'on'); ylabel(a(3),'spectral counts');
    
    % ------ PERCENT DIFFERENCE IN SPECTRUM CHANNELS
    plot(a(4),data.datenum(idx_ful_good),pdiff.spec_full(idx_ful_good,idx_wavelengths),'*','MarkerSize',3,'LineWidth',2);
    hold(a(4),'on'); grid(a(4),'on'); ylim(a(4),[-100 100])
    title(a(4),title_string),
    ylabel(a(4),'% Difference')
    
    %% Plot segments
    %plot(a(1),data_avg.datenum(idx_seg_good),data_avg.density(idx_seg_good),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(2),data_avg.datenum(idx_seg_good),data_avg.NO3_uM(idx_seg_good), '*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(3),data_avg.datenum(idx_seg_good),data_avg.spec_channels(idx_seg_good,idx_wavelengths),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(4),data_avg.datenum(idx_seg_good),pdiff.spec_seg(idx_seg_good,idx_wavelengths),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    
    a(1).XLim = xlimits;
    a(2).XLim = xlimits;
    a(3).XLim = xlimits;
    a(4).XLim = xlimits;
    a(2).XTickLabel = [];
    a(3).XTickLabel = [];
    a(4).XTickLabel = [];
    % link axes so when zoom all axes will change limits
    linkaxes([a(1) a(2) a(3) a(4)],'x')
  end %% function a = plot_segment_percent_diff 

%% function a = plot_segment_percent_diff_abs
  function a = plot_segment_percent_diff_abs
    title_string = 'Absorbance at 254nm and 350nm';
    % plot density, NO3, spectrum channel, percent differences and return axes
    a = makefig_subplots(1,4);% spectra
    xlimits = [data.datenum(1)-1 data.datenum(end)+1];
    idx_ful_good = find(data.dark ~= 0 & data.flag <= meta_proc.flag.not_evaluated); % & all(isfinite(data.spectrum_channels),2));
    idx_seg_good = find(data_avg.flag <= meta_proc.flag.not_evaluated);              % & all(isfinite(data_avg.spec_channels),2));
    %% Plot full resolution data
    % ------ DENSITY
    plot(a(1),data.datenum,data.density,'k-*','MarkerSize',2); % plot(a(1),data.datenum(idx_ful_good),data.density(idx_ful_good),'k-*','MarkerSize',2);
    hold(a(1),'on'); grid(a(1),'on'); 
    a(1).YLim = round(prctile(data.density(idx_ful_good),[0.02 99.99])); %ylim(a(1),[1023 1027])
    ylabel(a(1),'Density'); datetick(a(1),'Keepticks');
    
    % ------ RAW NO3 [MICROMOLAR]
    plot(a(2),data.datenum(idx_ful_good),data.NO3_uM(idx_ful_good),'k*','MarkerSize',3,'LineWidth',2);
    hold(a(2),'on'); grid(a(2),'on');ylabel(a(2),'NO_3 [\muM]');
    a(2).YLim = round(prctile(data.NO3_uM(idx_ful_good),[0.02 99.99])); %ylim(a(2),no3_limits) 
    
    % ------ ABSORBANCE AT 254 AND 350 
    plot(a(3),data.datenum(idx_ful_good),data.abs_254nm(idx_ful_good),'b*','MarkerSize',3,'LineWidth',2);
    hold(a(3),'on'); grid(a(3),'on'); ylabel(a(3),'Absorbance');
    plot(a(3),data.datenum(idx_ful_good),data.abs_350nm(idx_ful_good),'r*','MarkerSize',3,'LineWidth',2);
    
    % ------ PERCENT DIFFERENCE IN SPECTRUM CHANNELS
    plot(a(4),data.datenum(idx_ful_good),pdiff.abs_254nm_full(idx_ful_good),'b*','MarkerSize',3,'LineWidth',2);
    hold(a(4),'on'); grid(a(4),'on'); ylim(a(4),[-100 100])
    plot(a(4),data.datenum(idx_ful_good),pdiff.abs_350nm_full(idx_ful_good),'r*','MarkerSize',3,'LineWidth',2);
    title(a(4),title_string),
    ylabel(a(4),'% Difference')
    
    %% Plot segments
    %plot(a(1),data_avg.datenum(idx_seg_good),data_avg.density(idx_seg_good),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(2),data_avg.datenum(idx_seg_good),data_avg.NO3_uM(idx_seg_good), '*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(3),data_avg.datenum(idx_seg_good),data_avg.abs_254nm(idx_seg_good),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(3),data_avg.datenum(idx_seg_good),data_avg.abs_350nm(idx_seg_good),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(4),data_avg.datenum(idx_seg_good),pdiff.abs_254nm_seg(idx_seg_good),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);
    plot(a(4),data_avg.datenum(idx_seg_good),pdiff.abs_350nm_seg(idx_seg_good),'*','MarkerSize',2,'LineWidth',2,'Color',segment_clr);

    a(1).XLim = xlimits;
    a(2).XLim = xlimits;
    a(3).XLim = xlimits;
    a(4).XLim = xlimits;
    a(2).XTickLabel = [];
    a(3).XTickLabel = [];
    a(4).XTickLabel = [];
    % link axes so when zoom all axes will change limits
    linkaxes([a(1) a(2) a(3) a(4)],'x')
  end %% function a = plot_segment_percent_diff_abs
return
end