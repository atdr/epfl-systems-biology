function formatFig(w,h)
fig = gcf;
fig.PaperOrientation = 'landscape';
fig.PaperSize = [w h];
fig.PaperPosition = [0 0 w h];
fig.Renderer = 'Painters'; % for 3D plots
end