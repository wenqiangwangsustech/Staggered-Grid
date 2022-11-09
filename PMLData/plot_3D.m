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
soureceX =  150;
soureceY =  150;
soureceZ = 150;
dt = 0.0005;
filename = 'CPU StressXX';

    for j = start : step : Nt%20 :460
            fileName = [ 'stressXXSliceXX_XoY_', num2str( j ), '.txt' ];
            
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

% load curveLine1.txt;
% %hold on;
% plot( curveLine1 );
% title( 'damping = 0' )


%     fileName = [ 'snapshot_', num2str( 0 ),'_', num2str( 500 ), '.txt' ];
%    
%     data = load(fileName);
%     V = reshape( data, [ xLength, yLength, zLength ] );
%     [x,y,z] = meshgrid( 1 :xLength, 1 : yLength, 1 : zLength);
%     xslice = soureceX; yslice = soureceY; zslice = soureceZ;
%     %h = slice( x, y, z, V,xslice,yslice,zslice);
%     figure( 1 )
%       h1 = slice( x, y, z, V, xslice,[],[]); 
%       shading interp;
%     set(h1,'edgecolor','none');
%       figure( 2  )
%       h2 = slice( x, y, z, V, [],yslice,[]);
%       shading interp;
%     set(h2,'edgecolor','none');
%       figure( 3  )
%       h3 = slice( x, y, z, V, [],[],zslice);
%     %colormap( 'gray' )
%     shading interp;
%     set(h3,'edgecolor','none');
%     drawnow
% axis image
% %caxis(  [ -1e-7, 1e10-7] )
% colorbar






% for i = 100 : 100 : Nt
%     fileName = [ 'snapshot_0_', num2str(i ), '.txt' ];
%     data = load(fileName);
%     V = reshape( data, [ xLength, yLength, zLength ] );
%     [x,y,z] = meshgrid( 1 :xLength, 1 : yLength, 1 : zLength);
%     xslice = soureceX; yslice = soureceY; zslice = soureceZ;
%     h = slice( x, y, z, V,xslice,yslice,zslice);
%     %colormap( 'gray' )
%     shading interp;
%     set(h,'edgecolor','none');
%     drawnow
% end
% clear















% clear;
% close all;
% clc
% Nt = 300;
% xLength = 81;
% yLength = 81;
% zLength = 81;
% 
% % soureceX = xLength / 2;
% % soureceY =  yLength / 2;
% % soureceZ =  zLength / 2;
%  soureceX =  41;
% soureceY =  41;
% soureceZ =  41;
% figure( 1 );
%  dampingNameX =  'dampingX.txt' ;
%  dampingNameY =  'dampingY.txt' ;
%   dampingNameZ =  'dampingZ.txt' ;
% dataX = load(dampingNameX);
% dataY = load(dampingNameY);
% dataZ = load(dampingNameZ);
% V1 = reshape( dataX, [ xLength, yLength, zLength ] );
% V2 = reshape( dataY, [ xLength, yLength, zLength ] );
% V3 = reshape( dataZ, [ xLength, yLength, zLength ] );
% V = V1 + V2 + V3;
% [x,y,z] = meshgrid( 1 :xLength, 1 : yLength, 1 : zLength);
% xslice = soureceX; yslice = soureceY; zslice = soureceZ;
% h = slice( x, y, z, V,xslice,yslice,zslice);
% 
% %colormap( 'gray' )
% shading interp;
% set(h,'edgecolor','none');   
% %caxis(  [ -10, 10] ) 
% colorbar
% axis image
% drawnow







