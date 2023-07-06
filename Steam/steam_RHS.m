function dvdt = steam_RHS(~,z,par)
x=z(1);
v=z(2);
%precompute coefficients to visually declutter
%equation 9 will look something like:  A*vdot = Bv^2 + C* O(x)
sb_sc= par.sb / par.sc;
sb_sl= par.sb / par.sl;
A = par.xbar+x+ sb_sc*(par.ybar - sb_sc*x) + sb_sl*par.L ;
% fprintf("x= %f \n" , x);
assert( par.ybar - sb_sc*x> 0 ,"Liquid level in the conduit is negative");% assert that the liquid level in the conduit is positive
dyn_press = par.rho*par.g*(par.ybar-par.sb/par.sc*x -par.xbar-x);
% fprintf("x= %f \n P= %f \n" , x, dyn_press);
B = (-0.5)*(1-(sb_sc)^2);
C = -par.g * (par.xbar+x);
D = par.g * (par.ybar - (sb_sc*x));

vol = par.sb*(par.H-par.xbar-x);
assert(vol>0); % assert that the bubble trap volume is positive

E = (1/par.rho).*par.Fdudp(vol)*par.FdPdv(vol);
F = -par.alpha * (1- par.gamma)*(1 / par.rho) * par.Pa0;

dvdt =[v; (1./A)*(B*v^2 +C+D+E+F)];
end