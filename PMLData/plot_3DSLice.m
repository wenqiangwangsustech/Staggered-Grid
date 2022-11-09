clear;
close all;
clc
Nt = 1000;
xLength = 300;
yLength = 300;
zLength = 300;

step = 50;
start =950;
soureceX = xLength / 2;
soureceY =  yLength / 2;
soureceZ =  zLength / 2;

colorRange =  [ -1e-10, 1e-10 ];

h = figure;
title([ 'StressXX Middle Slice' ],'FontSize', 16);
dt = 0.0005;
for i = start : step : Nt
    
display( i );
dataXoY = load( ['snapshotSliceXX_XOY_', num2str( i ),'.txt' ]);
dataXoZ = load( ['snapshotSliceXX_XoZ_', num2str( i ),'.txt' ]);
dataYoZ = load( ['snapshotSliceXX_YoZ_', num2str( i ),'.txt' ]);
% dataXoY = load( ['velocitySliceXX_XoY_', num2str( i ),'.txt' ]);
% dataXoZ = load( ['velocitySliceXX_XoZ_', num2str( i ),'.txt' ]);
% dataYoZ = load( ['velocitySliceXX_YoZ_', num2str( i ),'.txt' ]);
% dataXoY = load( ['surfSliceXX_XOY_', num2str( i ),'.txt' ]);
% dataXoZ = load( ['surfSliceXX_XoZ_', num2str( i ),'.txt' ]);
% dataYoZ = load( ['surfSliceXX_YoZ_', num2str( i ),'.txt' ]);
%colorRange =  [ min( dataXoY ) / 3, max( dataXoY ) / 3];
subplot( 2,2, 1 );
pcolor( reshape( dataXoY, xLength, yLength) );
title([ 'XoY' , ' t = ',num2str( i * dt),'s'],'FontName','Times New Roman','FontSize' , 10);
shading interp;
caxis( colorRange ) 
colorbar
axis image
camlight
drawnow;

subplot( 2, 2, 2 );
%colorRange =  [ min( dataXoZ ) / 3, max( dataXoZ ) / 3];
pcolor( reshape( dataXoZ, xLength, yLength) );
title([ 'XoZ' , ' t = ',num2str( i * dt),'s'],'FontName','Times New Roman','FontSize',10 );
shading interp;
caxis(  colorRange ) 
colorbar
axis image
camlight
drawnow;

subplot( 2, 2, 3 );
%colorRange =  [ min( dataYoZ ) / 3, max( dataYoZ ) / 3];
pcolor( reshape( dataYoZ, xLength, yLength) );
title([ 'YoZ' , ' t = ',num2str( i * dt),'s'],'FontName','Times New Roman','FontSize',10 );
shading interp;
caxis( colorRange ) 
colorbar
axis image
camlight
drawnow;
drawGif( i, dt, 'Stress XX',  start, h )
end




