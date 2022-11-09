clear;
close all;
clc
Nt = 1500;
xLength = 100;
yLength = 200;
zLength = 300;

step = 100;
start = 500;
soureceX = xLength / 2;
soureceY =  yLength / 2;
soureceZ =  zLength / 2;

colorRange =  [ -1e-10, 1e-10 ];

figure;


for i = start : step : Nt
    
display( i );
% dataXoY = load( ['snapshotSliceXX_XOY_', num2str( i ),'.txt' ]);
% dataXoZ = load( ['snapshotSliceXX_XoZ_', num2str( i ),'.txt' ]);
% dataYoZ = load( ['snapshotSliceXX_YoZ_', num2str( i ),'.txt' ]);
% dataXoY = load( ['velocitySliceXX_XoY_', num2str( i ),'.txt' ]);
% dataXoZ = load( ['velocitySliceXX_XoZ_', num2str( i ),'.txt' ]);
% dataYoZ = load( ['velocitySliceXX_YoZ_', num2str( i ),'.txt' ]);
dataXoY = load( ['surfSliceXX_XOY_', num2str( i ),'.txt' ]);
dataXoZ = load( ['surfSliceXX_XoZ_', num2str( i ),'.txt' ]);
dataYoZ = load( ['surfSliceXX_YoZ_', num2str( i ),'.txt' ]);
%colorRange =  [ min( dataXoY ) / 3, max( dataXoY ) / 3];
subplot( 2,2, 1 );
pcolor( reshape( dataXoY, 300, 200 ) );
title( 'XoY' );
shading interp;
caxis( colorRange ) 
colorbar
axis image
camlight
drawnow;

subplot( 2, 2, 2 );
%colorRange =  [ min( dataXoZ ) / 3, max( dataXoZ ) / 3];
pcolor( reshape( dataXoZ, 300, 100 ) );
title( 'XoZ' );
shading interp;
caxis(  colorRange ) 
colorbar
axis image
camlight
drawnow;

subplot( 2, 2, 3 );
%colorRange =  [ min( dataYoZ ) / 3, max( dataYoZ ) / 3];
pcolor( reshape( dataYoZ, 200, 100 ) );
title( 'YoZ' );
shading interp;
caxis( colorRange ) 
colorbar
axis image
camlight
drawnow;

end




