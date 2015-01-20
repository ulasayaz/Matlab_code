function convertToEPS(imgName,type)

pic = imread(imgName,type);
[y x c] = size(pic)

figure('Units','Pixels','Resize','off',...
   'Position',[100 100 x y],'PaperUnit','points',...
   'PaperPosition',[0 0 x y]);
axes('position',[0 0 1 1]);
image(pic);
axis off

saveas(gcf,imgName,'epsc');