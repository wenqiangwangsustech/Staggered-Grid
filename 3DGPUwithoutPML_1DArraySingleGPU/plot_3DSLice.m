clear;
close all;
clc
Nt = 1000;
xLength = 500;
yLength = 500;
zLength = 500;

step = 30;

soureceX = xLength / 2;
soureceY =  yLength / 2;
soureceZ =  zLength / 2;

colorRange =  [ -0.001, 0.001];

figure;


for i = 10 : step : Nt
    
display( i );
dataXoY = load( ['stressXXSliceXX_XoY_', num2str( i ),'.txt' ]);
dataXoZ = load( ['stressXXSliceXX_XoZ_', num2str( i ),'.txt' ]);
dataYoZ = load( ['stressXXSliceXX_YoZ_', num2str( i ),'.txt' ]);
% dataXoY = load( ['velocityXSliceXX_XoY_', num2str( i ),'.txt' ]);
% dataXoZ = load( ['velocityXSliceXX_XoZ_', num2str( i ),'.txt' ]);
% dataYoZ = load( ['velocityXSliceXX_YoZ_', num2str( i ),'.txt' ]);

subplot( 2,2, 1 );
pcolor( reshape( dataXoY, xLength, yLength ) );
title( 'XoY' );
shading interp;
caxis( colorRange ) 
colorbar
axis image
camlight
drawnow;

subplot( 2, 2, 2 );
pcolor( reshape( dataXoZ, xLength, yLength ) );
title( 'XoZ' );
shading interp;
caxis(  colorRange ) 
colorbar
axis image
camlight
drawnow;

subplot( 2, 2, 3 );
pcolor( reshape( dataYoZ, xLength, yLength ) );
title( 'YoZ' );
shading interp;
caxis( colorRange ) 
colorbar
axis image
camlight
drawnow;

end




