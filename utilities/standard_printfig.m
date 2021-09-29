function standard_printfig(figname)
%% function standard_printfig(figname)
% print figure in standard way
% author: Brita Irving <bkirving@alaska.edu>
%%
fprintf(' saving figure to %s\n',figname)
image_type       = '-djpeg';
image_resolution = '-r250';
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(fig, 'InvertHardcopy', 'off')
print(figname,image_type,image_resolution)
end %% FUNCTION STANDARD_PRINTFIG(FIGNAME)