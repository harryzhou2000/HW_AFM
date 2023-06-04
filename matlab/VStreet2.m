syms a b k x y
assume([a,b],'positive');
% a = 1
% b = 1

Uk = 1/((b+y) * (1 + (k*a - x)^2/(b+y)^2));
Vk = simplify((k*a -x)/(b + y) * Uk);

Ukm = -y/(y^2 + (k*a - x)^2);
Vkm = -simplify((k * a - x)/(y^2 + (k*a - x)^2));

Uk_x = subs(diff(Uk,x),[x,y],[0,0])
Uk_y = subs(diff(Uk,y),[x,y],[0,0])

Vk_x = subs(diff(Vk,x),[x,y],[0,0])
Vk_y = subs(diff(Vk,y),[x,y],[0,0])


Ukm_x = subs(diff(Ukm,x),[x,y],[0,0])
Ukm_y = subs(diff(Ukm,y),[x,y],[0,0])

Vkm_x = subs(diff(Vkm,x),[x,y],[0,0])
Vkm_y = subs(diff(Vkm,y),[x,y],[0,0])


U_x = symsum(subs(Uk_x,k,(1 + 2*k)/2),k,0,inf) ...
     +symsum(subs(Uk_x,k,(1 + 2*k)/2),k,-inf,-1) 
 
U_y = symsum(subs(Uk_y,k,(1 + 2*k)/2),k,0,inf) ...
     +symsum(subs(Uk_y,k,(1 + 2*k)/2),k,-inf,-1) 
 
V_x = symsum(subs(Vk_x,k,(1 + 2*k)/2),k,0,inf) ...
     +symsum(subs(Vk_x,k,(1 + 2*k)/2),k,-inf,-1) 
 
V_y = symsum(subs(Vk_y,k,(1 + 2*k)/2),k,0,inf) ...
     +symsum(subs(Vk_y,k,(1 + 2*k)/2),k,-inf,-1) 
 
Um_x = symsum(Ukm_x,k,1,inf) ...
     +symsum(Ukm_x,k,-inf,-1) 
 
Um_y = symsum(Ukm_y,k,1,inf) ...
     +symsum(Ukm_y,k,-inf,-1) 
 
 
Vm_x = symsum(Vkm_x,k,1,inf) ...
     +symsum(Vkm_x,k,-inf,-1) 
 
Vm_y = symsum(Vkm_y,k,1,inf) ...
     +symsum(Vkm_y,k,-inf,-1) 
%%
UU_y = symfun(U_y+Um_y, [a,b])
VV_x = symfun(V_x+Vm_x, [a,b])
%%
double(U_y + Um_y)
double(V_x + Vm_x)
%%
vpa(VV_x(1,0.365))



