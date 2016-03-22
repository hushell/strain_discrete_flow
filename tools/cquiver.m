function cquiver(px,py)

assert(isequal(size(px),size(py)));

[ny,nx] = size(px);

x = 1:nx;
y = 1:ny;

mg = sqrt(px.^2+py.^2);

hold on
pcolor(x,y,mg);
%axis([-3 3 -3 3]);
axis([1 nx 1 ny]);
colormap((jet+white)/2);
colorbar
shading interp
quiver(x,y,px,py,2,'k');
axis off
hold off