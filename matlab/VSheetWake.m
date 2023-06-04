
gB = 1;
gR = 0.5;
xi0 = 1;
N = 360;
theta0s = linspace(0, 2*pi -1/N, N);

g0 = (gR + gB * sin(theta0s)) * 2 * pi * xi0 / N;
G0 = -2*pi * xi0 * gR;

xg0 = cos(theta0s) * xi0;
yg0 = sin(theta0s) * xi0;
xG0 = 0;
yG0 = 0;

gG0 = [g0,G0];
xgG0 = [xg0, xG0];
ygG0 = [yg0, yG0];





tEnd = 4;
dt = 0.0001;
tc = 0;
u = reshape([xgG0;ygG0],[],1);
%%
vout = VideoWriter('pVorts3.mp4', 'MPEG-4');
vout.FrameRate = 30;
vout.open();

for iter =1:1000000
    tnext= tc + dt;
    ifEnd = false;
    dtc = dt;
    if(tnext >= tEnd)
        dtc = tnext - tc;
        ifEnd = true;
        tnext = tEnd;
    end
    
    dudt0 = odeF(u,gG0);
    u1 = u + dudt0 * dt;
    dudt1 = odeF(u1,gG0);
    unew = u + dudt0 * dt /2  + dudt1 * dt/2;
    
    
    tc = tnext;
    u = unew;
    
    if(mod(iter,100) == 0 || ifEnd)
        figure(1); clf; cla;
        set(gcf,'Position',[100,100,600,600]);
        plotPs(u,gG0);
        axis equal;
        title(sprintf('t = %g, \\chi = %g', tc, gR));
        ylim([-2,2]);
        set(gca,'Position',[0.1300 0.1100 0.7750 0.8150]);
        set(gca,'OuterPosition',[0 0 1 1] );
        drawnow;
        fprintf("iter = %d\n",iter);
        
        
        frame = getframe(gcf);
        vout.writeVideo(frame);
    end
    if(ifEnd)
        break;
    end
end
vout.close();
%%
tmax = 1;
options = odeset('AbsTol',1e-4,'RelTol',1e-4,'InitialStep',0.01);
sol = ode113(@(t,v) odeF(v,gG0), [0,tmax], [xgG0;ygG0],options);


figure(1);clf; cla;
plotUana(sol.y(:,end), tmax, gR/gB);
plotPs( sol.y(:,end),gG0)
axis equal;
ylim([-2,2]);
title(sprintf('t = %g, \\chi = %g', tmax, gR));

L = legend('Analytic','PointVort');
%%
figure(2);
plotHs( sol.y(:,end), gG0)

%%



function [vxs,vys] = velo_pointVerts(xs,ys,gs)
dxs = xs - xs';
dys = ys - ys';
rs = (dxs.^2 + dys.^2);
rs = rs + eye(size(rs));
vxsa = (-dys./rs) .* (gs');
vysa = ( dxs./rs) .* (gs');
vxs = sum(vxsa,1)/2/pi;
vys = sum(vysa,1)/2/pi;
end

function vs = odeF(xs,gs)
xs = reshape( xs,2,[]);
[vxs,vys] = velo_pointVerts(xs(1,:), xs(2,:), gs);
vs = reshape([vxs;vys],[],1);
end

function plotPs(u,gGs)
xs = reshape(u,2,[]);
xc = xs(1,1:end-1);
yc = xs(2,1:end-1);
Xc = xs(1,end);
Yc = xs(2,end);
gs = gGs(1:end-1);

hold on;
plot(xc,yc, '-','DisplayName','PointVort');
mGs = max(abs(gs),[],'all');
scatter(xc,yc, 1, abs(gs))

hhmap = hot(256);
colormap(hhmap(end:-1:1,:));
c = colorbar;
c.Label.String = 'Relative Abs Vort Strength';
c.Position = [0.86,0.1100 0.04 0.8150];
caxis([-0.1 * mGs,inf]);

plot(Xc,Yc,'o');
hold off;


end

function plotHs(u,gGs)

xs = reshape(u,2,[]);
xc = xs(1,1:end-1);
yc = xs(2,1:end-1);
Xc = xs(1,end);
Yc = xs(2,end);
gs = gGs(1:end-1);

% histogram2(xc,yc);
[N,Xedges,Yedges] = histcounts2(xc,yc,[100,200]);
Xm = 0.5*(Xedges(2:end) + Xedges(1:end-1));
Ym = 0.5*(Yedges(2:end) + Yedges(1:end-1));
p = pcolor(Xm,Ym, N'/sum(N,'all'));
p.LineStyle = 'none';
axis xy;
colormap copper;

hold on;
plot(Xc,Yc,'o');
hold off;

end


function xi = xiiAna(theta,t, chi)

xi = 1 - 1/8 * t^2 * cos(2*theta) +...
    t^3 * (1/24*chi * sin(2*theta) + 1/32*cos(3*theta))...
    + t^4 * (1/96 *chi^2*cos(2*theta) - 5/256*chi*sin(3*theta) + 5/768*cos(2*theta) - 7/768*cos(4*theta));

end

function plotUana(u, t, chi)
    xs = reshape(u,2,[]);
Xc = xs(1,end);
Yc = xs(2,end);

thetas = linspace(0,2*pi, 1000);
xis = xiiAna(thetas,t,chi);
xs = cos(thetas) .* xis + Xc;
ys = sin(thetas) .* xis + Yc;

hold on;
plot(xs,ys,'DisplayName' ,'Analytic','LineWidth',2);



end



