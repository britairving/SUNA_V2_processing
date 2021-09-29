function standard_printfig_highrespng(figname)
%% function fig = makefig
% print figure in standard way
% author: Brita Irving <bkirving@alaska.edu>
%%
fprintf(' saving figure to %s\n',figname)
image_type       = '-dpng';
image_resolution = '-r350';
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(fig, 'InvertHardcopy', 'off')
print(figname,image_type,image_resolution)
end