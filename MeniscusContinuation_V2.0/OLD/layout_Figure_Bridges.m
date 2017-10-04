
%close(10);
gcf = figure(10);

set(gcf,'PaperUnits','centimeters')
xSize = 9.5; ySize = 8.5;
xLeft = 0; yTop = 0;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[1 1 xSize*60 ySize*60])

width = 3/xSize; height = 1.3/ySize;
x = .5/xSize; y = 1-(.3+1.3)/ySize; 
axes('position',[x y width height])
plot([0 0],[0,1.3],'-.k')
axis([0 3 0 1.3])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

width = 3/xSize; height = 2/ySize;
x = .5/xSize; y =1-(.3+1.3+.3+2)/ySize; 
axes('position',[x y width height])
plot([0 0],[0,2],'-.k')
axis([0 3 0 2])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

width = 3/xSize; height = 4/ySize;
x = .5/xSize; y =1-(.3+1.3+.3+2+.3+4)/ySize; 
axes('position',[x y width height])
plot([0 0],[0,4],'-.k')
axis([0 3 0 4])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

width = 4.5/xSize; height = 8/ySize;
x = (.5+3+.5)/xSize; y =.25/ySize; 
axes('position',[x y width height])
plot([0 0],[0,8],'-.k')
axis([0 4.5 0 8])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';




%x2 = x + width + 1/xSize;
%axes('position',[x y width height])