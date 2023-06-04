syms a b k
assume([a,b],'positive');

Uk = 1/(b * (1 + k^2*a^2/b^2));
Vk = k*a/b * Uk;



Uk_k = diff(Uk,k)
Uk_b = diff(Uk,b)
Vk_k = diff(Vk,k)
Vk_b = diff(Vk,b)
%%
U_k = symsum(subs(Uk_k,k,(1+2*k)/2),k, 0, inf) + ...
      symsum(subs(Uk_k,k,(1+2*k)/2),k, -inf, -1)

  
U_b = symsum(subs(Uk_b,k,(1+2*k)/2),k, 0, inf) + ...
      symsum(subs(Uk_b,k,(1+2*k)/2),k, -inf, -1)
  
V_k = symsum(subs(Vk_k,k,(1+2*k)/2),k, 0, inf) + ...
      symsum(subs(Vk_k,k,(1+2*k)/2),k, -inf, -1)
  
V_b = symsum(subs(Vk_b,k,(1+2*k)/2),k, 0, inf) + ...
      symsum(subs(Vk_b,k,(1+2*k)/2),k, -inf, -1)
  
V_kF = symfun(V_k, [a, b])

U_bF = symfun(U_b, [a, b])

%%
re = [];
Rs = linspace(0.01,5,100);
for R = Rs
 re(end+1) = double(U_bF(R,1) * V_kF(R,1));
end
plot(Rs,re)