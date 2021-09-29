function ctd = cnv2mat_mooring(ctd_file)
% CNV2MAT Reads single SeaBird ASCII .CNV file
%
%  Usage:   ctd=cnv2mat(ctd_dir);
%
%     Input:  ctd_file = path to cnv file
%
%     Output: lon = longitude in decimal degrees, West negative
%             lat = latitude in decimal degrees, North positive
%           gtime = Gregorian time vector in UTC
%            data = matrix containing all the columns of data in the .CNV file
%           names = string matrix containing the names and units of the columns
%         sensors = string matrix containing the names of the sensors
%
%  NOTE: How lon,lat and time are written to the header of the .CNV
%        file may vary with CTD setup.  For our .CNV files collected on
%        the Oceanus, the lat, lon & time info look like this:
%
%    * System UpLoad Time = Mar 30 1998 18:48:42
%    * NMEA Latitude = 42 32.15 N
%    * NMEA Longitude = 069 28.69 W
%    * NMEA UTC (Time) = 23:50:36
%
%  Modify the lat,lon and date string handling if your .CNV files are different.

%  4-8-98  Rich Signell (rsignell@usgs.gov)
%     incorporates ideas from code by Derek Fong & Peter Brickley
%  11-2018 Brita Irving (bkirving@alaska.edu)
%     update to read all cnv files

%% ctd into ctd directory and get list of all cnv files
% cd(ctd_dir)
% Open the .cnv file as read-only text
dbstop if error
%% build structure to store data
ctd = struct();
cnv_file = ctd_file;
cnv_file = char(cnv_file);
fid=fopen(cnv_file,'rt');
if fid == -1
  fprintf('could not read file %s\n',cnv_file)
  keyboard
else
  fprintf('Reading cnv file %s\n',cnv_file)
end
ctd.filename = cnv_file;
% Read the header.
% Start reading header lines of .CNV file,
% Stop at line that starts with '*END*'
%
% Pull out NMEA lat & lon along the way and look
% at the '# name' fields to see how many variables we have.
%
str='';
line_num = 0;
while ~contains(str,'*END*')
  str=fgetl(fid);
  if str == -1
    break
  else
    line_num = line_num + 1;
  end
  %-----------------------------------
  %
  %    Read the NMEA latitude string.  This may vary with CTD setup.
  %
  try
    
    if contains(str,'NMEA Lat')
      is=findstr(str,'=');
      isub=is+1:length(str);
      dm=sscanf(str(isub),'%f',2);
      if(findstr(str(isub),'N'))
        ctd.lat = dm(1)+dm(2)/60;
      else
        ctd.lat = -(dm(1)+dm(2)/60);
      end
      %-------------------------------
      %
      %    Read the NMEA longitude string.  This may vary with CTD setup.
      %
    elseif contains(str,'NMEA Lon')
      is=findstr(str,'=');
      isub=is+1:length(str);
      dm=sscanf(str(isub),'%f',2);
      if(findstr(str(isub),'E'));
        ctd.lon = dm(1)+dm(2)/60;
      else
        ctd.lon = -(dm(1)+dm(2)/60);
      end
      %------------------------
      %
      %    Read the 'System upload time' to get the date.
      %           This may vary with CTD setup.
      %
      %    I'm reading this in to get the date, since the NMEA time string
      %    does not contain date.  Unfortunately, the system upload time is
      %    in local time (here, EST), so I need to convert to UTC by adding
      %    5 hours (5/24 days).
      %
    elseif contains(str,'Station')
      is=strsplit(str,':');
      if contains(is{2},'Station')
        station = strrep(is{2},'Station','');
      else
        station = is{2};
      end
      ctd.station = strtrim(station);
    elseif contains(str,'System UpLoad')
      is=findstr(str,'=');
      %    pick apart date string and reassemble in DATEFORM type 0 form
      datstr=[str(is+6:is+7) '-' str(is+2:is+4) '-' str(is+9:is+12)];
      datstr=[datstr ' ' str(is+14:is+21)];
      %    convert datstr to Julian time, add 5 hours to convert from EST to GMT
      n=datenum(datstr)+5/24;
      gtime=datevec(n);
      %----------------------------
      %
      %    Read the NMEA TIME string.  This may vary with CTD setup.
      %
      %      replace the System upload time with the NMEA time
    elseif contains(str,'NMEA UTC')
      
      is=findstr(str,':');
      s = strsplit(str,'=');
      s2 = s{2};
      ctd.gtime = datenum(s2,'mmm dd yyyy HH:MM:SS');
      %------------------------------
      %
      %    Read the variable names & units into a cell array
      %
      
    elseif contains(str,'name')
      var=sscanf(str(7:10),'%d',1);
      var=var+1;  % .CNV file counts from 0, Matlab counts from 1
      %      stuff variable names into cell array
      names{var}=str;
      %------------------------------
      %
      %    Read the sensor spans
      %
      %     elseif contains(str,'span')
      %       depth_loc = find(contains(names,'Depth'));
      %       str2 = strsplit(names{depth_loc},'=');
      %       str2 = str2{1};
      %       str2 = strrep(str2,'name','span');
      %       if contains(str,str2)
      %         str3 = strsplit(str,',');
      %         ctd.botdepth = str2double(str3{2});
      %       end
      %       %------------------------------
      %       %
      %       %    Read the sensor names into a cell array
      %       %
    elseif contains(str,'# sensor')
      sens=sscanf(str(10:11),'%d',1);
      sens=sens+1;  % .CNV file counts from 0, Matlab counts from 1
      %      stuff sensor names into cell array
      sensors{sens}=str;
      %
      %  pick up bad flag value
    elseif contains(str,'# bad_flag')
      isub=13:length(str);
      bad_flag=sscanf(str(isub),'%g',1);
    end
  catch
    fprintf('stopped here\n')
    keyboard
  end
end
%==============================================
%
%  Done reading header.  Now read the data!
nvars=var;  %number of variables
% Read the data into one big matrix

data=fscanf(fid,'%f',[nvars inf]);
fclose(fid);

%
% Flag bad values with nan
%
ind=find(data==bad_flag);
data(ind)=data(ind)*nan;
%
% Flip data around so that each variable is a column
data=data.';
if size(data,1) > 1
  % Convert cell arrays of names to character matrices
  names=char(names);
  ctd.names_str = names;
  ctd.names     = {};
  ctd.units     = {};
  for nl = 1:size(names,1)
    lines = strsplit(names(nl,:),'=');
    line  = strsplit(lines{2},':');
    ctd.names(nl) = line(1);
    ctd.units(nl) = line(2);
  end
else
  try
    %fprintf('only read in a single like of data... see whats going on\n')
    % Try to read exact variable names
    fid=fopen(cnv_file,'rt');
    str = fgetl(fid);
    for nline = 1:line_num - 2
      str = fgetl(fid);
    end
    fclose(fid);
    str(1) = [];
    vars = strsplit(str,' ');
    if isempty(vars{1})
      vars(1) = [];
    end
    if any(contains(vars,'ParSbeox0Mm/Kg'))
      split_here = find(contains(vars,'ParSbeox0Mm/Kg'));
      vars1 = vars(1:split_here-1);
      vars2 = vars(split_here+1:end);
      vars = [vars1 'Par' 'Sbeox0Mm/Kg' vars2];
      vars_orig = vars;
      vars = strrep(vars,'/','_');
      vars = strrep(vars,'-','_');
      vars = strrep(vars,'é','e');
      vars = strrep(vars,':','_');
    end
    opt = detectImportOptions(cnv_file,'FileType','text','NumHeaderLines',line_num,'Delimiter',{' ' '\t'});
    opt.DataLine = line_num + 1;
    %   opt.LeadingDelimitersRule = 'ignore';
    %   opt.Whitespace = ' ';
    %   opt.Encoding = 'UTF-8';
    opt.ConsecutiveDelimitersRule = 'join';
    opt.LeadingDelimitersRule     = 'ignore';
    opt.VariableNamesLine = line_num-1;
    opt.VariableNames = vars;
    opt.VariableTypes(:) = {'double'};
    if any(contains(vars,'MMM')) % month in string format... e.g. AUG
      opt.VariableTypes(contains(vars,'MMM')) = {'char'};
    end
    if any(contains(vars,'HH_MM_SS')) % month in string format... e.g. AUG
      opt.VariableTypes(contains(vars,'HH_MM_SS')) = {'char'};
    end
    data = readtable(cnv_file,opt);
    
    % Update names & units
    %names=char(names);
    ctd.names_str = vars_orig;
    ctd.names     = vars;
    ctd.units     = {};
    
    for nl = 1:numel(vars)
      unit_line = char(names(contains(names,vars_orig{nl},'IgnoreCase',true)));
      try
        lines = strsplit(unit_line,'=');
        line  = strsplit(lines{2},':');
        %ctd.names(nl) = line(1);
        ctd.units(nl) = line(2);
      catch
        if any(strcmp(vars_orig{nl},{'DD' 'MMM' 'YYYY' 'HH:MM:SS'}))
          ctd.units(nl) = vars_orig(nl);
        else
          keyboard
        end
      end
    end % Loop through variables to pull out units
    try
      dates = [num2str(data.DD) repmat('-',size(data.YYYY)) char(data.MMM) repmat('-',size(data.YYYY)) num2str(data.YYYY) repmat(' ',size(data.YYYY)) char(data.HH_MM_SS)];
      data.datetime = datetime(dates,'InputFormat','dd-MMM-yyyy HH:mm:SS');
      data.datenum  = datenum(data.datetime);
      ctd.names = [ctd.names 'datetime' 'datenum'];
      ctd.units = [ctd.units 'dd-mmm-yyyy HH:MM:SS' 'days since January 0, 0000' ];
    catch
      fprintf('datenum and datetime not available..??\n')
      keyboard
    end
  catch
    fprintf('CTD file format not as expected, see whats going on...\n')
    keyboard
  end
end

% Save to ctd structure
ctd.data = data;
try
  % replace with standard variable names...
  for nl = 1:numel(ctd.names)
    varname0 = strtrim(ctd.names{nl});
    switch varname0
      case {'timeJV2' 'timeJ'};         varname = 'julday';
      case {'depSM' 'DepSM'};           varname = 'depth';
      case {'prDM' 'prdM' 'PrDM'};      varname = 'pressure';
      case {'tv290C' 'Tv290C'};         varname = 'temp';
      case {'sal00' 'Sal00'};           varname = 'salinity';
      case {'sigma-t00' 'Sigma_e00'};   varname = 'sigmat';
      case {'c0S/m' 'C0S_m'};           varname = 'cond';
      case {'potemp090C' 'Potemp090C'}; varname = 'pottemp';
      case 'flag';                      varname = 'flag';
      case 'CStarAt0';                  varname = 'beam_attenuation';
      case 'turbWETbb0';                varname = 'turbidity';
      case {'wetCDOM' 'WetCDOM'};       varname = 'cdom';
      case {'flECO-AFL' 'FlECO_AFL'};   varname = 'fluor_mgm3';
      case {'par' 'Par'};               varname = 'PAR';
      case {'V0' 'v0'};                 varname = 'voltage0';
      case {'V1' 'v1'};                 varname = 'voltage1';
      case 'sbeox0PS';                  varname = 'ox_sat';
      case {'oxsolMm/Kg' 'OxsatMm_Kg'}; varname = 'ox_sat_umol_kg';
      case {'sbox0Mm/Kg' 'sbeox0Mm/Kg' 'Sbeox0Mm_Kg'};  varname = 'ox_umol_kg';
      otherwise
        varname = strrep(varname0,'-','_');
        varname = strrep(varname,'/','_');
        varname = strrep(varname,'é','');
    end
    ctd.names{nl} = varname;
    ctd.(varname) = ctd.data(:,nl);
    if istable(ctd.(varname))
      ctd.(varname) = table2array(ctd.(varname));
    end
    % prDM: Pressure, Digiquartz [db]
    % t090C: Temperature [ITS-90, deg C]
    % t190C: Temperature, 2 [ITS-90, deg C]
    % c0S/m: Conductivity [S/m]
    % c1S/m: Conductivity, 2 [S/m]
    % v0: Voltage 0
    % v1: Voltage 1
    % v2: Voltage 2
    % v3: Voltage 3
    % v4: Voltage 4
    % v5: Voltage 5
    % v6: Voltage 6
    % v7: Voltage 7
    % flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]
    % CStarAt0: Beam Attenuation, WET Labs C-Star [1/m]
    % CStarTr0: Beam Transmission, WET Labs C-Star [%]
    % par: PAR/Irradiance, Biospherical/Licor
    % sbox0Mm/Kg: Oxygen, SBE 43 [umol/kg]
    % sbox1Mm/Kg: Oxygen, SBE 43, 2 [umol/kg]
    % sbeox0Mm/Kg: Oxygen, SBE 43 [umol/kg]
    % sbeox1Mm/Kg: Oxygen, SBE 43, 2 [umol/kg]
    % sbeox0PS: Oxygen, SBE 43 [% saturation]
    % sbeox1PS: Oxygen, SBE 43, 2 [% saturation]
    % altM: Altimeter [m]
    % timeJ: Julian Days
    % sal00: Salinity, Practical [PSU]
    % sal11: Salinity, Practical, 2 [PSU]
    % sigma-t00: Density [sigma-t, kg/m^3 ]
    % sigma-t11: Density, 2 [sigma-t, kg/m^3 ]
    % upoly0: Upoly 0, Deep-SUNA
    % upoly1: Upoly 1, LISST
  end
  
  if exist('gtime','var') && ismember('julday',ctd.names)
    juldays = ctd.data(:,ismember(ctd.names,'julday'));
    if max(juldays > 365)
      ctd.datenum = datenum(gtime(1)-1,0,juldays);
    else
      ctd.datenum = datenum(gtime(1),0,juldays);
    end
    ctd.datetime = datetime(ctd.datenum,'ConvertFrom','datenum');
    ctd.names = [ctd.names 'datenum' 'datetime'];
    ctd.units = [ctd.units 'days since January 0, 0000' 'dd-mmm-yyyy HH:MM:SS'];
  end
  
  
  
  
  %   isfield(ctd,'
  
  %   % calculate the density using sw_dens.m
  %   % dens = sw_dens(S,T,P)
  %   P = data(:,1);
  %   T = data(:,4);
  %   S = data(:,10);
  %
  %   dens = sw_dens(S,T,P);
  
  %   % Plot Data
  %   % simplify file name
  %   st1 = strsplit(cnv_file, '.');
  %   st2 = strrep(st1(1),'_','-');
  %   test = char(st2);
  %   strtitle =[test ':' ' ' num2str(lat) char(176) 'N' ' ' num2str(lon) char(176) 'E'];
  %   % only plot to 1000 m / 1000 db
  %   if max(P) > 1000
  %     pmax = P <= 1000;
  %     plot(T(pmax),-P(pmax))
  %   end
  %   if max(P) < 1000
  %     plot(T,-P)
  %   end
  %   ylabel(' Pressure [db]')
  %   xlabel(' Temperature [C]')
  %   title(strtitle)
  %   saveas(gcf,cnv_file,'jpg')
  
catch
  keyboard
end
