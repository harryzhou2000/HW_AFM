sigma = 10;
b = 8/3;

fLorentz = @(t,y,sigma, b, r) [sigma*(-y(1) + y(2)); r*y(1)-y(2)-y(1)*y(3); -b*y(3) + y(1)*y(2)];

%% chaos
rng(1);
rcur = 25;
nsamp = 5;
range = [-1,1] * 10;
if(rcur < 1)
    base = [0;0;0];
    yinits = base + rand(3,nsamp) *(range(2) - range(1)) + range(1);
else
    base1 = [sqrt(b*(rcur-1)), sqrt(b*(rcur-1)), rcur - 1]';
    base2 = [-sqrt(b*(rcur-1)), -sqrt(b*(rcur-1)), rcur - 1]';
    yinits = [base1 + rand(3,nsamp) *(range(2) - range(1)) + range(1),...
        base2 + rand(3,nsamp) *(range(2) - range(1)) + range(1)];
end


figure(1);cla; clf; hold on; ax1 = gca;
figure(2);cla; clf; hold on; ax2 = gca;

for is = 1:size(yinits,2)
    
    sol = fsolveLorentz(fLorentz, yinits(:,is), 50, sigma, b, rcur);
    
    xends = sol.y(1,sol.x >0.5);
    
%     fmt = "r";
%     if all(xends > 0) || all (xends < 0)
%         fmt = "b";
%     end
%     plotLorentz (ax1, sol, fmt);
%     plotLorentzX(ax2, sol, fmt);
    
    plotLorentz (ax1, sol);
    plotLorentzX(ax2, sol);
    
end

figure(1);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
view(3);

figure(2);
xlabel('\tau');
ylabel('X');
grid on;


%% limit cycle
rng(1);
rcur = 20;
nsamp = 20;
range = [-8, 8];
if(rcur < 1)
    base = [0;0;0];
    yinits = base + rand(3,nsamp) *(range(2) - range(1)) + range(1);
else
    base1 = [sqrt(b*(rcur-1)), sqrt(b*(rcur-1)), rcur - 1]';
    base2 = [-sqrt(b*(rcur-1)), -sqrt(b*(rcur-1)), rcur - 1]';
    yinits = [base1 + rand(3,nsamp) *(range(2) - range(1)) + range(1),...
        base2 + rand(3,nsamp) *(range(2) - range(1)) + range(1)];
end


figure(1);cla; clf; hold on; ax1 = gca;
figure(2);cla; clf; hold on; ax2 = gca;

for is = 1:size(yinits,2)
    
    sol = fsolveLorentz(fLorentz, yinits(:,is), 50, sigma, b, rcur);
    
    xends = sol.y(1,sol.x >0.5);
    
    fmt = "r";
    if all(xends > 0) || all (xends < 0)
        fmt = "b";
    end
    plotLorentz (ax1, sol, fmt);
    plotLorentzX(ax2, sol, fmt);
    
    % plotLorentz (ax1, sol);
    % plotLorentzX(ax2, sol);
    
end

figure(1);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
view(3);

figure(2);
xlabel('\tau');
ylabel('X');
grid on;


%% simple equil
rng(1);
rcur = 2;
nsamp = 8;
range = [-1, 1];
if(rcur < 1)
    base = [0;0;0];
    yinits = base + rand(3,nsamp) *(range(2) - range(1)) + range(1);
else
    base1 = [sqrt(b*(rcur-1)), sqrt(b*(rcur-1)), rcur - 1]';
    base2 = [-sqrt(b*(rcur-1)), -sqrt(b*(rcur-1)), rcur - 1]';
    yinits = [base1 + rand(3,nsamp) *(range(2) - range(1)) + range(1),...
        base2 + rand(3,nsamp) *(range(2) - range(1)) + range(1)];
end


figure(1);cla; clf; hold on; ax1 = gca;
figure(2);cla; clf; hold on; ax2 = gca;

for is = 1:size(yinits,2)
    
    sol = fsolveLorentz(fLorentz, yinits(:,is), 50, sigma, b, rcur);
    
    xends = sol.y(1,sol.x >0.5);
    
%     fmt = "r";
%     if all(xends > 0) || all (xends < 0)
%         fmt = "b";
%     end
%     plotLorentz (ax1, sol, fmt);
%     plotLorentzX(ax2, sol, fmt);
    
    plotLorentz (ax1, sol);
    plotLorentzX(ax2, sol);
    
end

figure(1);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
view(3);

figure(2);
xlabel('\tau');
ylabel('X');
grid on;


%%
figure(1);
print(gcf, sprintf("XYZ_r_%g.png", rcur), '-dpng', '-r300');
figure(2);
print(gcf, sprintf("XT_r_%g.png", rcur), '-dpng', '-r300');



%%
function sol = fsolveLorentz(f, yinit, tmax,sigma, b, r)

opts = odeset('AbsTol',1e-12, 'RelTol', 1e-12);

sol = ode113(@(t,y) f(t,y,sigma,b,r), [0,tmax], yinit, opts);





end

function p = plotLorentz(a, sol, fmt)
if( nargin == 2)
    p = plot3(a, sol.y(1,:), sol.y(2,:), sol.y(3,:));
else
    p = plot3(a, sol.y(1,:), sol.y(2,:), sol.y(3,:), fmt);
end

end
function p = plotLorentzX(a, sol, fmt)
if( nargin == 2)
    p = plot(a, sol.x, sol.y(1,:));
else
    p = plot(a, sol.x, sol.y(1,:), fmt);
end
end
