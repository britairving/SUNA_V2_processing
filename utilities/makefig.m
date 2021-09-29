function fig = makefig
%% function fig = makefig
%
% Description:
%   Make a new figure. Gives user ability to customize figure font, size,
%   color, etc. 
%   If there are two monitors, tries to use the second screen.
%
%
% author: Brita Irving <bkirving@alaska.edu>
 %% Create standardized figure (on second monitor, if possible)
try 
  FontChoice = 'Times';
  axiscolor  = [0.9 0.9 0.9];
  % Create new figure
  n = 1;
  while ishandle(n)
    n = n + 1;
  end
  fig = figure(n);
  % First, set units to normalized so can set screen size
  fig.Units = 'normalized';
  % change size to just a bit smaller than screen size
  fig.OuterPosition = [0.02 0.04 0.90 0.90];
  % Next, move figure to second monitor if possible
  % If the get(0, 'MonitorPositions') function worked, and there are 2
  % monitors, position the figure so it is on the second screen.
  try
    % Get monitor positions in pixel size
    MP = get(0, 'MonitorPositions');
    if size(MP,1) > 1
      % Change units to Pixels because that is what get(0,[]) function returns
      fig.Units = 'pixels';
      % Now change position to second monitor
      fig.Position(1) = MP(2,1)+5;
      pause(0.2); % for some reason.. need to pause to change screen
      if ~isequal(fig.Position(1),MP(2,1)+5)
        fig.Position(1) = MP(2,1)+5;
      end
    end
  catch
    % plot will stay in first screen
  end

  fig.Units = 'Inches';
  fig.PaperPositionMode = 'auto';
  fig.PaperOrientation = 'portrait';
  % axes defaults
  set(fig,'DefaultTextFontname',    FontChoice);
  set(fig,'defaultTextFontSize',    16);
  set(fig,'DefaultTextFontName',    FontChoice);
  set(fig,'DefaultAxesFontName',    FontChoice);
  set(fig,'DefaultAxesFontSize',    20);
  set(fig,'DefaultAxesFontWeight',  'Bold');
  set(fig, 'Color', 'w');
  set(fig,'DefaultAxesLineStyleOrder','-|--|-.|:');
  
  tomLnCmap = [ 0 0 1; 1 0 0; 1 0.6 0; 0.85 0 0.95; 0 .95 0.95; 0.8 0.8 0; 0.8 0.5 0; 0.5 0 0; 0 0 0.5; 0 0.65 0.65];
  %               b;     r;   orange;      ~m;        ~bc;          ~y;    dorange;      dr;     db;       dc;
  % no greens as Tom is R-G color blind.
  set(fig,'DefaultAxesColorOrder',tomLnCmap);
  % set axis color
  ax = gca;
  % ydir = reverse is the default because want to show ocean depth from 0
  % (top) to full depth (bottom)
  set(ax,'YDir','reverse','XAxisLocation','bottom','Layer','top','Color',axiscolor);

catch % Create a simple figure using Matlab's default
  % This will happen if any of the above fails
  fig = figure;
  keyboard
end % try/catch
end % MAIN FUNCTION FIG = MAKEFIG