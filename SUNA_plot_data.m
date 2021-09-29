function SUNA_plot_data(data_proc,meta_proc,cfg)
%% function SUNA_plot_data(data_proc,meta_proc,data_avg,cfg)
%
%
%
%
cal = meta_proc.NO3_uM_TCSS.caldata;
%% Plot: Spectrogram
cmap = cmocean('thermal');
spec_limits = [0 60000];
time_limits = [round(data_proc.datenum(1)-1) round(data_proc.datenum(end)+1)];

for nn = 1:2
  makefig;
  a1 = subplot(3,1,1:2);
  a2 = subplot(3,1,3);
  if nn == 1
    % show all data
    idx_good     = find(data_proc.dark ~= 0);
    title_string = [strrep(cfg.project,'_','\_') ' | SUNA ' cfg.mooring.SN.SUNA ' | Spectra & Nitrate before automatic QC'];
    ylims = round(prctile(data_proc.NO3_uM(idx_good),[0.02 99.99])); %
    if ylims(1) > 0; ylims(1) = 0; end
    if max(data_proc.NO3_uM_cor(idx_good)) > ylims(2) 
      ylims(2) = ceil(max(data_proc.NO3_uM_cor(idx_good)));
    end
    figname = fullfile(cfg.datadir,[cfg.project '_spectra_beforeQC']);
  else
    % Only show data with flags less than or equal to not_evaluated
    idx_good     = find(data_proc.dark ~= 0 & data_proc.flag <= meta_proc.flag.not_evaluated);
    title_string = [strrep(cfg.project,'_','\_') ' | SUNA ' cfg.mooring.SN.SUNA ' | Spectra & Nitrate after automatic QC'];
    figname = fullfile(cfg.datadir,[cfg.project '_spectra_afterQC']);
    %	ylims = [0 ylims(2)];
  end
  % Plot spectrum without dark counts
  fprintf('\nimagesc: '); tic
  % hp = pcolor(a1,data_proc.datenum(idx_good),cal.Wavelength,data_proc.spectrum_channels(idx_good,:)');
  hp = imagesc(a1,data_proc.datenum(idx_good),cal.Wavelength,data_proc.spectrum_channels(idx_good,:)');
  shading(a1,'flat');
  a1.YLabel.String = 'Wavelength [nm]';
  a1.YDir = 'normal';
  cb = colorbar(a1);
  cb.Limits = spec_limits;
  a1.CLim   = cb.Limits;
  a1.XLim   = time_limits;
  datetick(a1,'x','yyyy-mm','keeplimits');

  colormap(a1,cmap);
  a1.XTickLabel = [];
  toc
  
  % basic plot to show corrected values
  plot(a2,data_proc.datenum(idx_good),data_proc.NO3_uM(idx_good),'k-','MarkerSize',2,'LineWidth',1,'DisplayName','NO_3 Raw')
  hold(a2,'on'); grid(a2,'on');
  plot(a2,data_proc.datenum(idx_good),data_proc.NO3_uM_cor(idx_good),'b-','MarkerSize',2,'LineWidth',2,'DisplayName','NO_3 Corrected')
  % legend(a2,'show');
  ylabel(a2,'NO_3 [\muM]');
  a2.XLim   = time_limits;
  datetick(a2,'x','yyyy-mm','keeplimits');
  hl = legend(a2,'show'); hl.FontSize = 12;
  a2.YLim = ylims;

  a2.Position(3) = a1.Position(3);
  pause(4);
  a2.Position(3) = a1.Position(3);
  % align all axes so can easily spot drop outs
  a1.Title.String  = title_string;
  linkaxes([a1 a2],'x');
  
  
 if cfg.save_figures
   standard_printfig_highrespng(figname);
 end

end