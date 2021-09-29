function cals = calculate_co2sys_ph(type,cfg,cals,mtch)
%% FUNCTION CALCULATE_CO2SYS_PH
%
%  Syntax:
%    cals = calculate_co2sys_ph(type,cfg,cals,mtch)
%
%  Inputs:
%    type | 'calc' or 'insitu' 
%    cfg  | structure containing details about processing
%    cals | table containing calibration cast data
%    mtch | required for type='insitu'. Structure with insitu T,S,P
%
%  Description:
%    Use CO2SYS toolbox to calculate pH from available parameters, or
%    caluclate pH at insitu conditions. 
%
%  References:
%    https://github.com/jamesorr/CO2SYS-MATLAB
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Default coefficients for CO2sys
if isfield(cfg,'CO2SYS')
  phscale = cfg.CO2SYS.pHScale;
  k1k2    = cfg.CO2SYS.K1K2;
  kso4    = cfg.CO2SYS.KSO4;
else
  phscale = 1;  %1 = Total scale
  k1k2    = 10; %10 = Lueker et al, 2000	T:    2-35  S: 19-43. Total scale. Real seawater.
  kso4    = 1;  %1 = KSO4 of Dickson 1990a   & TB of Uppstrom 1974  (PREFERRED)
end

%% initialize new pH & pCO2 variables
switch type
  case 'calc' 
    cals.pH_calc   = nan(size(cals.pH));
    cals.pCO2_calc = nan(size(cals.pH));
  case 'insitu'
    cals.pH_insitu          = nan(size(cals.pH));
    cals.pCO2_insitu        = nan(size(cals.pH));
    cals.pH_insitu_pressure = nan(size(cals.pH));
end

% Loop through each sample
for nsamp = 1:size(cals,1)
  cast = cals(nsamp,:);
  %fprintf('Calculating in-situ pH using CO2SYS and nearby discrete samples | %s %s %d\n',char(cast.CRUISE), char(cast.StationID),cast.Niskin)
  ph_insitu = 1;
  if isfinite(cast.pH) && strcmp(type,'calc')
    % Discrete sample available, do not need to calculate
    continue
  else
    % Pull out silicate if unavailable
    if isfinite(cast.Sil); sil = cast.Sil;
    elseif isfinite(cast.Sil_uM); sil = round(cast.Sil_uM .* 1000 ./cast.density,3); % [umol/L] * 1000[L/m3] * 1/[kg/m3]
    else; sil = nanmean(data_avg.SIest);
    end
    % Pull out phosphate if unavailable
    if isfinite(cast.PO4); phos = cast.PO4;
    elseif isfinite(cast.PO4_uM); phos = round(cast.PO4_uM .* 1000 ./cast.density,3); % [umol/L] * 1000[L/m3] * 1/[kg/m3]
    else; phos = nanmean(data_avg.PO4est);
    end
    % Pull out phosphate if unavailable
    if isfinite(cast.TA); TA = cast.TA;        TA_u = cast.TA_u;
    elseif isfinite(cast.alk);  TA = cast.alk; %TA_u = cast.alk_u;
    else; TA = NaN;
    end
    % Pull out DIC if unavailable
    if isfinite(cast.DIC); DIC = cast.DIC;        DIC_u = cast.DIC_u;
    elseif isfinite(cast.TCO2);  DIC = cast.TCO2;% DIC_u = cast.TCO2_u;
    else; DIC = NaN;
    end
    if isequal(DIC,-999); DIC = NaN; end
    % Pull out pH if unavailable
    if isfinite(cast.pH); pH = cast.pH; pH_u = cast.pH_u;
    elseif isfinite(cast.pH_calc); pH = cast.pH_calc;
    else; pH = NaN;
    end
    
    % Calculate pH for discrete sample
    if strcmp(type,'calc') && isfinite(DIC) && isfinite(TA)
      [DATA_cal]=CO2SYS(DIC,TA,2,1,cast.Salinity,cast.Temp, cast.Temp,...
        cast.Pressure,cast.Pressure,sil,phos, phscale,k1k2,kso4);
      % Pull out pH and pCO2
      cals.pH_calc(nsamp)   = round(DATA_cal(18),4);
      cals.pCO2_calc(nsamp) = round(DATA_cal(19),4);
    % Calculate pH at insitu conditions
    elseif strcmp(type,'insitu')
      cast_name = strrep(cast.StationID,'.','_');
      cast_name = char(strrep(cast_name,'-','_'));
      if isfield(mtch,cast_name) && ~isempty(mtch.(cast_name).idx_ctd)
        % Calculate pH for discrete sample
        ph = [];
        if isfinite(DIC) && isfinite(TA) && isfinite(pH)
          [DATA_cal1]=CO2SYS(DIC,TA,2,1,cast.Salinity,cast.Temp, mtch.(cast_name).temperature,...
            cast.Pressure,mtch.(cast_name).pressure,sil,phos, phscale,k1k2,kso4);
          [DATA_cal2]=CO2SYS(pH,TA,3,1,cast.Salinity,cast.Temp, mtch.(cast_name).temperature,...
            cast.Pressure,mtch.(cast_name).pressure,sil,phos, phscale,k1k2,kso4);
          [DATA_cal3]=CO2SYS(DIC,pH,2,3,cast.Salinity,cast.Temp, mtch.(cast_name).temperature,...
            cast.Pressure,mtch.(cast_name).pressure,sil,phos, phscale,k1k2,kso4);
          ph   = round(nanmean([DATA_cal1(18) DATA_cal2(18) DATA_cal3(18)]),4);
          pco2 = round(nanmean([DATA_cal1(19) DATA_cal2(19) DATA_cal3(19)]),4);
        elseif isfinite(DIC) && isfinite(pH)
          [DATA_cal]=CO2SYS(DIC,pH,2,3,cast.Salinity,cast.Temp, mtch.(cast_name).temperature,...
            cast.Pressure,mtch.(cast_name).pressure,sil,phos, phscale,k1k2,kso4);
          ph   = round(DATA_cal(18),4);
          pco2 = round(DATA_cal(19),4);
        elseif isfinite(TA) && isfinite(pH)
          [DATA_cal]=CO2SYS(pH,TA,3,1,cast.Salinity,cast.Temp, mtch.(cast_name).temperature,...
            cast.Pressure,mtch.(cast_name).pressure,sil,phos, phscale,k1k2,kso4);
          ph   = round(DATA_cal(18),4);
          pco2 = round(DATA_cal(19),4);
        else
          ph_insitu = 0;
          pco2      = 0;
        end
        if ph_insitu
          cals.pH_insitu(nsamp)   = ph;
          cals.pCO2_insitu(nsamp) = pco2;
          cals.pH_insitu_pressure(nsamp) = mtch.(cast_name).pressure;
          cals.pH_insitu_salinity(nsamp) = mtch.(cast_name).salinity;
          cals.pH_insitu_temperature(nsamp) = mtch.(cast_name).temperature;
        end
      end
    end
  end
end %% Loop through calibration samples and calculate pH
end %% FUNCTION CALCULATE_CO2SYS_PH

