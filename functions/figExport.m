function figExport(w,h,name)
formatFig(w,h)
print(gcf, '-dpdf', [pwd '/figures/' name '.pdf']);
end