function cfg = mooring_config(cfg)
%FUNCTION mooring_config 
%
%  Syntax:
%    cfg = mooring_config(cfg)
%  
%  Description:
%    Defines global metadata such as creator, publisher, contributor, etc.
%    Calls individual mooring configuration script based on cfg.project
%    that defines basic metadata such as deployment date, recovery date,
%    sensor serial numbers, data processing flags, and loads relevant
%    calibration cast data.
%
%  References: 
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%% Define basic metadata
cfg.metadata.global.featureType         = 'timeSeries';
cfg.metadata.global.platform            = 'mooring';
cfg.metadata.global.platform_vocabulary = 'http://mmisw.org/ont/ioos/platform';
cfg.metadata.global.asset_type          = 'network';      % network = A network of stations defined above from https://ioos.github.io/conventions-for-observing-asset-identifiers/ioos-assets-v1-0.html
cfg.metadata.global.cdm_data_type       = 'Station';
cfg.metadata.global.conventions         = 'CF-1.8, ACDD-1.3';

%% Insitution, creator, contributor, and publisher information
cfg.metadata.global.contributor_name      = 'Seth Danielson, Brita Irving, Peter Shipton';
cfg.metadata.global.contributor_role      = 'Principal Investigator, Data Manager, Mooring Technician';
cfg.metadata.global.institution           = 'University of Alaska Fairbanks';
cfg.metadata.global.creator_institution   = 'University of Alaska Fairbanks, College of Fisheries and Ocean Science';
cfg.metadata.global.creator_type          = 'group';
cfg.metadata.global.creator_address       = '2150 Koyukuk Drive';
cfg.metadata.global.creator_city          = 'Fairbanks';
cfg.metadata.global.creator_country       = 'USA';
cfg.metadata.global.creator_phone         = '907-474-7834';
cfg.metadata.global.creator_sector        = 'academic';
cfg.metadata.global.creator_state         = 'Alaska';
cfg.metadata.global.creator_zipcode       = '99775-7220';
% Publisher information
cfg.metadata.global.publisher_name        = 'Brita Irving';
cfg.metadata.global.publisher_institution = 'University of Alaska Fairbanks, International Arctic Research Center';
cfg.metadata.global.publisher_url         = 'github.com/britairving, britairving.com, claudinehauri.com';
cfg.metadata.global.publisher_orcid       = 'https://orcid.org/0000-0002-9474-6203';
cfg.metadata.global.publisher_address     = '2160 Koyukuk Drive';
cfg.metadata.global.publisher_city        = 'Fairbanks';
cfg.metadata.global.publisher_country     = 'USA';
cfg.metadata.global.publisher_phone       = '907-474-5966';
cfg.metadata.global.publisher_state       = 'Alaska';
cfg.metadata.global.publisher_zipcode     = '99775-7220';

%% Pull out mooring specific information
% ----------- Chukchi Sea Ecosystem Observatory ---------------------------
if contains(cfg.project,'CEO')
  cfg = CEO_config(cfg);
% ----------- Gulf of Alaska Ecosystem Observatory ------------------------
elseif contains(cfg.project,'GEO')
  cfg = GEO_config(cfg);
else
  error('Project name not recognized, need to configure')
end

%% Define data quality flags
% Code | Value | Definition
% 1    | Good | Passed documented QC tests
% 2    | Not evaluated, not available, or unknown | Used for data when no QC test was performed, or the information on quality is not available
% 3    | Questionable | Failed non-critical documented metric or subjective test
% 4    | Bad | Failed critical documented QC test(s) or as assigned by the data provider
% 5    | Estimate | Cell value were interpolated, extrapolated, or otherwise estimated
% 6    | Below detection limit | Value is below the detection limit of the analytical methods applied
% 7    | Outside instrument specification | Value is within the extended measuring range
% 9    | Missing Data | Used as placeholder when data are missing
cfg.flag = struct();
cfg.flag.good                    = 1;
cfg.flag.not_evaluated           = 2;
cfg.flag.questionable            = 3;
cfg.flag.bad                     = 4;
% cfg.flag.estimate                = 5;
% cfg.flag.below_detection_limit   = 6;
cfg.flag.outside_instrument_spec = 7;
cfg.flag.missing_data            = 9;

%% % Auto generate flag metadata -- these are variable attributes, NOT global
flags = fieldnames(cfg.flag);
cfg.metadata.flag_meanings = '';
cfg.metadata.flag_values   = [];
for n = 1:numel(flags)
  cfg.metadata.flag_meanings = [cfg.metadata.flag_meanings ' ' flags{n}];
  cfg.metadata.flag_values   = [cfg.metadata.flag_values, cfg.flag.(flags{n}) ];
end
% get rid of leading white space
cfg.metadata.flag_meanings = strtrim(cfg.metadata.flag_meanings);