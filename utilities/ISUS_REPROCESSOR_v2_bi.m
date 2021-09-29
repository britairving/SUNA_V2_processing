function [NTR_full,rmse_full] = ISUS_REPROCESSOR_v2_bi(data_proc,meta_proc,cal)
% Reprocesses absorption spectra sampled by the Satlantic ISUS/SUNA using
% co-located CTD data. 
%
%
% The centerpiece of the reprocessing, i.e. Sakamoto et al.'s algorithm
% (2009), is implemented in the file ISUS_REPROCESSOR_v2_2.m. You can also
% find an example of a ISUS/SUNA calibration file, and how to translate it
% into a .mat file that the reprocessing algorithm understands. The
% original inputs to the ISUS_REPROCESSOR_v2_2.m function were:
%     CHANNEL: the absorption data, one value per channel per sample,
%     CHANNEL_DF: the dark frames as measured in the same ISUS data file,
%     in-situ temperature T,
%     salinity S, and
%     a string calfile that supplies the name of the calibration .mat file. 
% 
% Adapted from: https://github.com/poplarShift/isus-suna-reprocessing
%   Randelhoff, A., Fer, I., Sundfjord, A., Tremblay, J.-É., & Reigstad, M.
%   (2016). Vertical Fluxes of Nitrate in the Seasonal Nitracline of the
%   Atlantic Sector of the Arctic Ocean. Journal of Geophysical Research:
%   Oceans, 121(7), 5282–5295. https://doi.org/10.1002/2016JC011779
%
% Equations based on Sakamoto et al., 2009 
%   Sakamoto, C. M., K. S. Johnson, and L. J. Coletti (2009), Improved
%   algorithm for the computation of nitrate concentrations in seawater
%   using an in situ ultraviolet spectrophotometer, Limnol. Oceanogr.
%   Methods, 7, 132–143.
%
% v2.2: added spec structure in output...
%
% achim@npolar.no
% Aug 2015
%% Notes on calibration files
% From Sea-Bird Tech Support 
% ________________________________________________________________________
% QUESTION: In the cal files is the 3rd column, titled "NO3", the
%  concentration, or absorption? Can you please elaborate on it's 
%  derivation and units?
% ANSWER: The 3rd column, titled "NO3" is the ratio of UV absorbance
%  at each wavelength with the reference nitrate standard specified in the
%  calibration range for each SUNA.
% ________________________________________________________________________
% QUESTION: I've looked at  various .cal files and other than "Reference",
%  the values never change (see attached .cal files), is that correct?
% ANSWER: If you're seeing multiple calibration files where only the
%  reference column changed, that means that these were calibrations where 
%  a DIW reference update was performed by the user.  When a full 
%  calibration is performed at Sea-bird, all of the columns will be updated.
% ________________________________________________________________________

%% Limit to only good data 
% Pull out good data indicies
data = data_proc;
if ~isempty(meta_proc) 
  % Ignore both flagged data and dark counts
  idx_good = find(data.flag <= meta_proc.flag.not_evaluated & data.dark ~= 0);
else
  % ONLY ignore dark counts
  idx_good = find(data.dark ~= 0);
end
% 
fields = fieldnames(data);
for nf = 1:numel(fields)
  data.(fields{nf}) = data.(fields{nf})(idx_good,:);
end

%% Just pass NTR and rmse as output arguments, do not pass out and spec
% Setting limit_output to 1 significantly increases speed!
limit_output = 1;
  
%% Define coefficients
% A,B,C,D coefficients for temperature dependency of ESW on sample
% temperature [Sakamoto et.al. 2009]
coef.A =  1.1500276 ;
coef.B =  0.02840;
coef.C = -0.3101349;
coef.D =  0.001222 ;

%% Pull out wavelengths between 217 and 240nm
idL = find(cal.Wavelength>217 & cal.Wavelength<240);
fit_wavelength = cal.Wavelength(idL);

% % wavelength
% LL=cal.lambda;
% position vector 217-240 nm

% % calibration temperature
% Tcal=20;
% % ref spectrum from DIW
% IR=cal.ref_DIW;
% % NO3 spectrum  
% % RE: Seabird "The 3rd column, titled "NO3" is the ratio of UV absorbance
% % at each wavelength with the reference nitrate standard specified in the
% % calibration range for each SUNA.
% ENO3 = cal.ENO3 ;


%% Initialize variables that don't change in loop
NTR  = nan(size(data.spectrum_channels,1),1);
rmse = nan(size(data.spectrum_channels,1),1);
if ~limit_output
  Apmean = nan(size(data.spectrum_channels,1),1);
end

% --- COEFFICIENT MATRIX
XX = [cal.NO3 cal.Wavelength] ;

% coeff matrix incl. constant offset
LXX = [ones(size(cal.Wavelength)) XX] ;

%% Loop through rows of data (1 row = 1 timestamp)
for k= 1:length(NTR) %119600:220000%
    try progressbar(k/length(NTR));end
    
    % --- CALCULATE QUANTITIES (SPECTRA):
    %   I:  intensity
    %   ID: dark intensity
    %   ASM: measured absorbance
    %   ESW: estimated seawater extinction (full waveband)
    %   ASE: estimated seawater absorption (full waveband)
    %   Ap:  absorption spectrum corrected for seawater absorbance (incl. temp, full waveband)
    
    % signal intensity
    I = data.spectrum_channels(k,:)';
    % dark current
    %ID = data.dark(k ,:)';
    ID = repmat(data.dark(k),size(I));
    % measured absorbance [Sakamoto et al., 2009 Equation 1]
    ASM = -log10((I-ID)./(cal.Reference-ID)) ;
    
    % in situ temperature
    Tis = data.temperature(k);
    %     if isnan(Tis)
    %         Tis=2;
    %     end
    % extinction coeff. of seawater at meas.temperature, extended beyond
    % Sakamoto et al's range... normalized to 1 PSU !
    % -- calfile ESW, scaled to in-situ temp
    % ESW  = cal.ESW(:) .* (coef.A+coef.B*Tis) ./ (coef.A+coef.B*T_CAL) .* exp (coef.D*(Tis-T_CAL) .* (cal.Wavelength-210)) ;
    % -- sakamoto's fit  [Sakamoto et al., 2009 Equation 4]
    ESW = (coef.A+coef.B*Tis) .* exp( (coef.C+coef.D*Tis).*(cal.Wavelength-210)) / 35;
    % expected absorption spectrum due to seawater  [Sakamoto et al., 2009 Equation 6]
    ASE = ESW * data.salinity(k) ;
    % absorption spectrum corrected for seawater absorbance (incl. temp)
    % [Sakamoto et al., 2009 Equation 7]
    Ap = ASM - ASE ;
   
    % coefficient matrix in fitting waveband (idL)
    X = XX(idL,:) ;
    NN = size(XX,2)+1 ;
      
    if sum(~isnan(Ap))>length(Ap)/2 % only fit if sufficiently many data points exist
        % linear regression
        p = robustfit(X,Ap(idL),'ols') ; % ordinary least squares...
        NTR_tmp = p(2);
        %NTR(k,1) = p(2);
        Afit = LXX * p; % fitted absorption
        rmse_tmp = sqrt(nanmean((Afit(idL)-Ap(idL)).^2)) ; % RMS error of fit
        %rmse(k,1) = sqrt(nanmean((Afit(idL)-Ap(idL)).^2)) ; % RMS error of fit
        if ~limit_output
          ANO3 = cal.NO3 * p(2) ; % NO3 absorption
          pout(k,:)=p;
          index_bg = setdiff(1:NN,2)' ;
          Abg = LXX(:,index_bg) * p(index_bg); % background absorption
        end
    else
      NTR_tmp  = NaN;
      rmse_tmp = NaN;
      %NTR(k,1)= NaN;
      if ~limit_output
        Afit = nan(size(Ap));
        Abg = Afit ;
        ANO3 = Afit ;
        Apmean(k,1)=NaN;
        p=nan(1,NN);
      end
    end
    %% September 2020 - do not store full and fit spec because takes too long!
    if ~limit_output
      % record some spectral info for later plotting
      spec.full.Ap  (k,:) = Ap  ;
      spec.full.ASE (k,:) = ASE ;
      spec.full.ASM (k,:) = ASM ;
      spec.full.ANO3(k,:) = ANO3;
      spec.full.Afit(k,:) = Afit;
      spec.full.Abg (k,:) = Abg ;
      spec.full.I   (k,:) = I   ;
      spec.full.ID  (k,:) = ID  ;
      
      spec.fit.Ap  (k,:) = Ap  (idL) ;
      spec.fit.ASE (k,:) = ASE (idL) ;
      spec.fit.ASM (k,:) = ASM (idL) ;
      spec.fit.ANO3(k,:) = ANO3(idL) ;
      spec.fit.Afit(k,:) = Afit(idL) ;
      spec.fit.Abg (k,:) = Abg (idL) ;
      spec.fit.I   (k,:) = I   (idL) ;
      spec.fit.ID  (k,:) = ID  (idL) ;
      spec.p (k,:) = p ; % linear regression, ordinary least squares...
    end
    % Do not store value if not real
    if isreal(NTR_tmp)
      NTR(k)  = NTR_tmp;
      rmse(k) = rmse_tmp;
    end
end

%% Reindex to full resolution
NTR_full  = nan(size(data_proc.NO3_uM));
rmse_full = nan(size(data_proc.NO3_uM));
NTR_full(idx_good)  = NTR;
rmse_full(idx_good) = rmse;

%% September 2020 just pass NTR and rmse as direct output arguments
if ~limit_output
  spec.full.cal.Reference (k,:) = cal.Reference      ;
  spec.fit.cal.Reference  (k,:) = cal.Reference(idL) ;
  
  spec.full.lambda = cal.Wavelength';
  spec.fit.lambda  = fit_wavelength';
  
  %create header
  %HEADER=create_readme(mfilename,'') ;
  %out.HEADER      = HEADER;
  out.HEADER      = mfilename;
  out.NTR         = NTR_full;
  out.temperature = data_proc.temperature;
  out.salinity    = data_proc.salinity;
  out.rmse        = rmse_full;
end

fprintf('\n')

end %% MAIN FUNCTION