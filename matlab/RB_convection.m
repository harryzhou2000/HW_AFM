a = 1.4;
psi0 = 1;
xs = linspace(0,2/a,501);
ys = linspace(0,1,501);
zs = [-1,1];


[xm,ym, zm] = meshgrid(xs,ys, zs);

um = -pi*psi0*sin(a*pi*xm).*cos(pi*ym);
vm = a*pi*psi0*cos(a*pi*xm).*sin(pi*ym);
sizem1 = size(xm,1);
sizem2 = size(xm,2);

startx = linspace(0,2/a,51);
starty = 0 * startx + 0.5;

%%

cla;
hold on;
% slines = stream2(xm,ym,um,vm,startx, starty);
% slobj = streamline(slines,'arrowsmode','arrows');
streamslice(xm,ym,zm, um,vm,vm*0,[],[],0)
axis equal;

patch('XData',[0,0,2/a,2/a],'YData',[0,1,1,0],'FaceColor','none')
xlabel('x');
ylabel('y');
% quiver(xm,ym,um,vm)