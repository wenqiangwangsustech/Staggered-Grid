clear all;
a = 10;
b = 10;
c = 10;

ax= 10;
ax= 10;
az= 10;

Epsilon = 0.1;   %Von-Karman

l=50;     %m
v0=3500;  
length = 300;
halfLength = length / 2;
Nx = -halfLength : 1 : halfLength;
Ny = -halfLength : 1 : halfLength;
Nz = -halfLength : 1 : halfLength;

kx = Nx / l;
ky = Ny / l;
kz = Nz / l;

for i = 1 : length
    for j = 1 : length
        for k = 1 : length
            fm( i, j, k ) = ( a * b * c ) / ( 1 + ( kx( i ) * kx( i ) * a * a )^ax + ( ky( j ) * ky( j ) * b * b )^az + ( kz( k ) * kz( k ) * b * b )^az );
            %c( i, j, k ) = ( a * b * c ) / fm( i, j, k );
%          c(i,j)=exp(-1.0*sqrt(x(j)*x(j)/a/a+y(i)*y(i)/b/b));
        end
    end
end

figure( 1 )
[x,y,z] = meshgrid( 1 : length, 1 : length, 1 : length);

xslice = halfLength; 
yslice = halfLength; 
zslice = halfLength;

h = slice( x, y, z, fm,xslice,yslice,zslice)
set(h,'edgecolor','none');

Theta = 2 * pi * rand( length, length, length );
Phi =sqrt( fm );

for i = 1 : length
    for j = 1 : length
        for k = 1 : length
           Psi( i, j, k ) = Phi( i, j, k ) * exp( Theta( i, j, k ) * sqrt( -1 ) );
        end
    end
end

figure( 2 )
[x,y,z] = meshgrid( 1 : length, 1 : length, 1 : length);
xslice = halfLength; yslice = halfLength; zslice = halfLength;
h = slice( x, y, z, abs( Psi ),xslice,yslice,zslice)
set(h,'edgecolor','none');

ImagRandom = ifftn( Psi );
RealRandom=abs( ImagRandom  );
figure( 3 )
[x,y,z] = meshgrid( 1 : length, 1 : length, 1 : length);
xslice = halfLength; yslice = halfLength; zslice = halfLength;
h = slice( x, y, z, RealRandom,xslice,yslice,zslice)
set(h,'edgecolor','none');

u = sum( RealRandom( : )  ) / ( length * length * length );
for i = 1 : length
    for j = 1 : length
        for k = 1 : length
           Delta0( i, j, k ) = RealRandom( i, j, k )  - u ;
           SquareDelta0( i, j, k ) = Delta0( i, j, k ) * Delta0( i, j, k );
        end
    end
end

d = sqrt( sum( SquareDelta0( : ) ) / ( length * length * length ) );
Delta = 0.1 * Delta0 / d;

v0 = 2000;
for i = 1 : length
    for j = 1 : length
        for k = 1 : length
           V( i, j, k ) = v0 * ( Delta( i, j, k ) + 1 );
        end
    end
end
figure( 4 )
[x,y,z] = meshgrid( 1 : length, 1 : length, 1 : length);
xslice = 1; yslice = 1; zslice = 1;
h = slice( x, y, z, V,xslice,yslice,zslice)
set(h,'edgecolor','none');
hold on;
xslice = length; yslice = length; zslice = length;
h = slice( x, y, z, V,xslice,yslice,zslice)
set(h,'edgecolor','none');

    
Fid=fopen('V.dat','wb');
cnt=fwrite(Fid,V,'float');
fclose(Fid);
