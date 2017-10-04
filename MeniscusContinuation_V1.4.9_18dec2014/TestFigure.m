
%close(10);
gcf = figure(10);

set(gcf,'PaperUnits','centimeters')
xSize = 14; ySize = 7.5;
xLeft = 0; yTop = 0;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[1 1 xSize*60 ySize*60])

width = 1.5/xSize; height = 3/ySize;
x = .35/xSize; y = 1-(.5+3)/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 1.5 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

width = 4/xSize; height = 3/ySize;
x = (.35+1.5+.5)/xSize; y =1-(.5+3)/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 4 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';

width = 7/xSize; height = 3/ySize;
x = (.35+1.5+.5+4+.5)/xSize; y =1-(.5+3)/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 7 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';


width = 13/xSize; height = 3/ySize;
x = .5/xSize; y =.5/ySize; 
axes('position',[x y width height])
plot([0 0],[-.5,2.5],'-.k')
axis([0 13 -.5 2.5])
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';





%x2 = x + width + 1/xSize;
%axes('position',[x y width height])