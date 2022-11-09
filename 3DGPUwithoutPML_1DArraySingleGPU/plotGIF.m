close all;
clear;
clc

Nt = 990;
dt = 0.0005;
nx = 500;
ny = 500;
start = 10;
step = 10;

fileName{ 1 } = 'stressXXSliceXX_XoY_';
fileName{ 2 } = 'stressXXSliceXX_XoZ_';
fileName{ 3 } = 'stressXXSliceXX_YoZ_';

fileName{ 4 } = 'stressXYSliceXX_XoY_';
fileName{ 5 } = 'stressXYSliceXX_XoZ_';
fileName{ 6 } = 'stressXYSliceXX_YoZ_';

fileName{ 7 } = 'velocityXSliceXX_XoY_';
fileName{ 8 } = 'velocityXSliceXX_XoZ_';
fileName{ 9 } = 'velocityXSliceXX_YoZ_';



filename{ 1 } = 'stressXX XoY';
filename{ 2 } = 'stressXX XoZ';
filename{ 3 } = 'stressXX YoZ';

filename{ 4 } = 'stressXY XoY';
filename{ 5 } = 'stressXY XoZ';
filename{ 6 } = 'stressXY YoZ';

filename{ 7 } = 'velocityX XoY';
filename{ 8 } = 'velocityX XoZ';
filename{ 9 } = 'velocityX YoZ';



color_range = [ -0.0001, 0.0001];
deltaT = 0.0005;

for n = 1 : 9
h = figure( n );
    for i = start : step : Nt
        data = load( [fileName{ n }, num2str( i ),'.txt' ]);
        drawGif( data, i,deltaT, filename{ n }, color_range, nx, ny, start, h );
    end
close all;
end