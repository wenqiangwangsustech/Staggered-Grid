clear;
close all;
clc
Nt = 1500;
xLength = 300;
yLength = 300;
zLength = 300;
start = 50;
step = 50;
% soureceX = xLength / 2;
% soureceY =  yLength / 2;
% soureceZ =  zLength / 2;
soureceX =  300;
soureceY =  300;
soureceZ =  300;
dt = 0.0005;
filename = 'CPU StressXX';

    for j = start : step : Nt%20 :460
            fileName = [ 'snapshot_0_', num2str( j ), '.txt' ];
            
            data = load(fileName);
            V = reshape( data, [ xLength, yLength, zLength ] );
            [x,y,z] = meshgrid( 1 :xLength, 1 : yLength, 1 : zLength);
            xslice = soureceX; yslice = soureceY; zslice = soureceZ;
           h = slice( x, y, z, V,xslice,yslice,zslice);
            %pcolor( squeeze( V( soureceX, :, : )  ));
            %colormap( 'gray' )
            shading interp;
            set(h,'edgecolor','none');   
           %caxis(  [ -1e-2, 1e-2] ) 
            colorbar
            camlight
            axis image
            drawnow
             drawGif3D( j, dt, filename,  start, fig )
            
           
    end








