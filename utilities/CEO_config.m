function cfg = CEO_config(cfg)
% FUNCTION CEO_CONFIG
%
%  Syntax:
%    cfg = CEO_config(cfg)
%  
%  Description:
%    Build initial metadata by outlining project details, what variables,
%    calibration files, processing steps, etc etc
%
%  See also:
%    SEDIT
%    EDIT
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Define paths - USER INPUT REQUIRED
if ispc
  cfg.path.chukchiobs = 'D:\MATLAB\hauri_lab\Chukchi_modeling\observations';
  cfg.path.calcasts   = 'D:\HauriLab\CEOmooring\CalibrationSamples\';
else
  cfg.path.chukchiobs = '/Users/bkirving/Documents/MATLAB/hauri_lab/Chukchi_modeling/observations';
  cfg.path.calcasts   = '/Volumes/CEO/HauriLab/CEOmooring/CalibrationSamples/';
end
cfg.path.calfile = fullfile(cfg.path.calcasts,'CEO_CalibrationSamples.xlsx');

%% Define mooring basics
cfg.mooring.name = 'CEO'; % short hand to use in scripts

%% Define global attributes
cfg.metadata.global.mooring_name = 'Chukchi Sea Ecosystem Observatory';
cfg.metadata.global.mooring_code = 'CEO2';
% Basic metadata
cfg.metadata.global.sea_name    = 'Chukchi Sea';
cfg.metadata.global.summary     = 'A multi-institutional consortium maintaining year-round Arctic marine ecosystem monitoring & process studies. The observatory is located on the southern flank of Hanna Shoal, 70 miles offshore, at 71.6 °N and 161.5 °W.';
cfg.metadata.global.geospatial_lat_min = 71.5;
cfg.metadata.global.geospatial_lat_max = 71.7; 
cfg.metadata.global.geospatial_lat_units = 'degrees_north';
cfg.metadata.global.geospatial_lon_min = -161.6;
cfg.metadata.global.geospatial_lon_max = -161.4;
cfg.metadata.global.geospatial_lon_units = 'degrees_east';
cfg.metadata.global.program                   = 'Chukchi Sea Ecosystem Observatory';
cfg.metadata.global.references                = 'http://research.cfos.uaf.edu/ceo/index.html, https://doi.org/10.5194/os-14-1423-2018';
cfg.metadata.global.ioos_regional_association = 'Alaska Ocean Observing System';
cfg.metadata.global.land_acknowledgment       = 'Data were collected on the ancestral waters of the Inupiat people. We recognize their unique relationship with this place, and are grateful for their stewardship past, present, and future. We strive to uphold Indigenous sovereignty and be in good relations with the original peoples and with this place.';
cfg.metadata.global.acknowledgment            = 'This project would not be possible without the availability of donated ship time for mooring turn-arounds each year, for which we especially thank chief scientists KatrinIken(2015 and 2017), RussHopcroft(2016), and Carin Ashjian (2018).';

%% Define deployment info and sensors
% Update this as necessary
% Instrument information can be found on 
%  1) AXIOM Research Workspace | 1426 Chukchi Ecosystem Monitoring
%     <https://researchworkspace.com/project/302408/folder/5225532/mooring_diagrams>
%  2) Google drive folder | 1426 Chukchi Ecosystem Monitoring
%     <https://drive.google.com/drive/folders/1WISYD64eqQL_10WFieH3nWmkEZx9bemr?usp=sharing>
switch cfg.year
  case 2014
    cfg.mooring.deploydate  = datetime(2014,09,20,22,48,00);
    cfg.mooring.recoverdate = datetime(2015,08,20,00,00,00);
    cfg.mooring.latitude    = 71.5997;   %71deg 35.980min N
    cfg.mooring.longitude   = -161.5054; %161deg 30.321min W
    cfg.mooring.SN.SBE16    = '4604';
    cfg.mooring.SN.triplet  = '516';
    cfg.mooring.SN.PAR      = '3059';
    cfg.mooring.SN.AZFP     = '55063';
  case 2015 
    cfg.mooring.deploydate  = datetime(2015,08,20,16,53,00);
    cfg.mooring.recoverdate = datetime(2016,08,08,01,01,01);
    cfg.mooring.latitude    = 71.59979;  %71deg 35.9874min N
    cfg.mooring.longitude   = -161.5261; %161deg 31.566min W
    cfg.mooring.SN.SBE16    = '6670';
    cfg.mooring.SN.triplet  = '';
    cfg.mooring.SN.PAR      = '';
    cfg.mooring.SN.SUNA     = '562';
    cfg.mooring.SN.AZFP     = '55082';
  case 2016
    cfg.mooring.deploydate  = datetime(2016,08,04,04,33,14);
    cfg.mooring.recoverdate = datetime(2017,08,14,21,30,00);
    cfg.mooring.latitude    = 71.5996;   %71deg 35.976min N
    cfg.mooring.longitude   = -161.5184; %161deg 31.621min W
    cfg.mooring.SN.SBE16    = '4604';
    cfg.mooring.SN.triplet  = '1417';
    cfg.mooring.SN.PAR      = '70612';
    cfg.mooring.SN.SUNA     = '801';
    cfg.mooring.SN.HydroC   = 'CO2-0416-001';
    cfg.mooring.SN.SeaFET   = '380';
    cfg.mooring.SN.SBE63    = '1340';
    cfg.mooring.SN.SBE37    = '14507';
    cfg.mooring.SN.AZFP     = '55063';
  case 2017
    cfg.mooring.deploydate  = datetime(2017,08,15,01,03,00);
    cfg.mooring.recoverdate = datetime(2018,08,05,18,03,14);
    cfg.mooring.latitude    = 71.5997;   %71deg 35.98min N
    cfg.mooring.longitude   = -161.5189; %161deg 31.65min W
    cfg.mooring.SN.SBE16    = '6670';
    cfg.mooring.SN.triplet  = '3956';
    cfg.mooring.SN.PAR      = '70592';
    cfg.mooring.SN.SUNA     = '562';
    cfg.mooring.SN.HydroC   = 'CO2-1216-001';
    cfg.mooring.SN.SeaFET   = '448';
    cfg.mooring.SN.SBE63    = '1578';
    cfg.mooring.SN.SBE37    = '15353';
    cfg.mooring.SN.AZFP     = '55082';
  case 2018
    cfg.mooring.deploydate  = datetime(2018,08,06,05,00,00);
    cfg.mooring.recoverdate = datetime(2019,08,20,21,54,57);
    cfg.mooring.latitude    = 71.5999;   %71deg 35.997min N
    cfg.mooring.longitude   = -161.5281; %161deg 31.683min W
    cfg.mooring.SN.SBE16    = '4604';
    cfg.mooring.SN.triplet  = '1417';
    cfg.mooring.SN.PAR      = '70612';
    cfg.mooring.SN.SUNA     = '801';
    cfg.mooring.SN.HydroC   = 'CO2-0416-001';
    cfg.mooring.SN.SeaFET   = '380';
    cfg.mooring.SN.SBE63    = '1340';
    cfg.mooring.SN.SBE37    = '14507';
    cfg.mooring.SN.AZFP     = '55063';
  case 2019
    cfg.mooring.deploydate  = datetime(2019,08,19,13,15,00);
    cfg.mooring.recoverdate = datetime(2020,10,05,21,00,00);
    cfg.mooring.latitude    = 71.59966;    %71deg 35.9796min N
    cfg.mooring.longitude   = -161.527463; %161deg 31.64778min W
    cfg.mooring.SN.SBE16    = '6670';
    cfg.mooring.SN.triplet  = '3956';
    cfg.mooring.SN.PAR      = '70592';
    cfg.mooring.SN.SBE43    = '3135';
    cfg.mooring.SN.SUNA     = '562';
    cfg.mooring.SN.HydroC   = 'CO2-1216-001';
    cfg.mooring.SN.AZFP     = '55082';
  case 2020
    cfg.mooring.deploydate  = datetime(2020,10,7,00,48,00);
    cfg.mooring.recoverdate = datetime(2021);
    cfg.mooring.latitude    = 71.5999;   %71deg 35.9947min N
    cfg.mooring.longitude   = -161.5259; %161deg 31.5553min W
    cfg.mooring.SN.SBE16    = '';
    cfg.mooring.SN.triplet  = '6099';
    cfg.mooring.SN.PAR      = '534';
    cfg.mooring.SN.SUNA     = '840';
    cfg.mooring.SN.HydroC   = 'CO2-0416-001';
    cfg.mooring.SN.SBE63    = '1578';
    cfg.mooring.SN.SBE37    = '21696';
    cfg.mooring.SN.AZFP     = '55063';
  otherwise
    error('year not recognized or set up yet')
end

%% Read calibration cast data
cfg = CEO_read_calibration_cast_data(cfg);

%% Update metadata with specifics defined for each project
cfg.metadata.global.time_coverage_start = datestr(cfg.mooring.deploydate);
cfg.metadata.global.time_coverage_end   = datestr(cfg.mooring.recoverdate);


end