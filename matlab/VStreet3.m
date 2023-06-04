syms a b k kappa z0
assume([a,b],'real')
assume([a,b],'positive')
assume(kappa, 'real')
% assume(kappa, 'integer')
assume(~in(1/2 - (b*1i)/a, 'integer') & ~in((b*1i)/a + 1/2, 'integer') &...
    in(a,'positive') & in(b,'positive'))

wUk = 1/1i  * 1/(-k*a);
wLk = 1/1i * 1/(-k*a + b*1i);

w = symsum(wUk + subs(wUk,k,-k),k,1,inf) ...
    -symsum(subs(wLk,k,(1+2*k)/2) + subs(wLk,k,(-1-2*k)/2),0,inf) ...
    
ewUk = 1/1i * (-1)/(-k)^2 * (...
    1 - exp(1i*kappa *(k )));
ewLk = 1/1i * (-1)/(-k + b*1i)^2 * (...
    1 - exp(1i*kappa *(k  - 1i*b)));


ew = symsum(ewUk + subs(ewUk,k,-k),k,1,inf) ...
    -symsum(subs(ewLk,k,(1+2*k)/2) + subs(ewLk,k,(-1-2*k)/2),0,inf)
    
%%
Few = @(aa,bb,kk) symsum(subs(ewUk + subs(ewUk,k,-k),[a,b,kappa],[aa,bb,kk]),k,1,inf) ...
    -symsum(subs(subs(ewLk,k,(1+2*k)/2) + subs(ewLk,k,(-1-2*k)/2),[a,b,kappa],[aa,bb,kk]),0,inf)
