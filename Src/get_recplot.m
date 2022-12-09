function pp=get_recplot(xc,yc,sigma,down,up,step)
figure
pp=scatter3(xc,yc,sigma,90,sigma,'filled');
view([0 90])
colormap jet
colorb = colorbar;
caxis([down up])
colorb.Ticks=down:step:up;
colorb.TickLabelInterpreter='latex';
colorb.Label.String = '$\delta\sigma$ $(S/m)$';
colorb.FontSize=14;
colorb.Label.FontSize=16;
colorb.Label.Interpreter = 'latex';
axis off
end