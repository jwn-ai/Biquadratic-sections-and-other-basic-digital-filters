clear all;

%(z-z1)(z-z1*)/(z^2+z+1) complex poles, real zeros
% b2 = 1;b1=0;b0=-1;
% a2 = 1;a1 = 1;a0=1;
% [r p k] =residue([b2 b1 b0],[a2 a1 a0])


syms  w0 rho r n theta real;
syms z;

p = r*exp(i*w0);
z1 = rho*exp(i*theta);
pc = conj(p);
z1c = conj(z1)

%B = ((p-z1)*(p-z1c))/((p-pc)*p);
B = ((p-z1)*(p-z1c))*((pc-p)*pc);
C = (p-pc)*p;
simplify(C*conj(C))
Bc = conj(B);

exp1 = B*(z-p)+Bc*(z-pc);
b = simplify(exp1,'Steps',100)

den = (z-p)*(z-pc);
a = simplify(den,'Steps',100)

H = simplify(b/a,'Steps',150)


H = simplify(H/(C*conj(C)))

h = simplify(expand(B*r^n*exp(i*w0*n) + Bc*r^n*exp(-i*w0*n)),'Steps', 500)


H0 = ( z*(1-rho^2/r^2) - 2*r*cos(w0) +2*rho*cos(theta))./


syms w real


s = exp(i*w);
sc = conj(s);

H0 = (s-z2)*(s-z1)*(sc-z1)*(sc-z2);


H = simplify(H0,'Steps',100)

Hd0 = (s-p)*(s-pc)*(sc-pc)*(sc-p);

Hd = simplify(Hd0,'Steps',150)
