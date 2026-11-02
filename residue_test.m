

%(z^2-1)/(z^2+z+1) complex poles, real zeros
b2 = 1;b1=0;b0=-1;
a2 = 1;a1 = 1;a0=1;
[r p k] =residue([b2 b1 b0],[a2 a1 a0])


syms  w0 z1 z2 r n real;
syms z B Bc;

p = r*exp(i*w0);

pc = conj(p);

B = (p-z1)*(p-z2)/((p-pc)*p);
Bc = conj(B);

exp1 = B*(z-p)+Bc*(z-pc);
b = simplify(exp1,'Steps',100)

den = (z-p)*(z-pc);
a = simplify(den,'Steps',100)

H = simplify(b/a,'Steps',150)
t = subs(H,z1=0)
subs(t,z2=r*cos(w))

B = (p-z1)*(p-z2)/((pc-p)*pc);
Bc = conj(B);

h = simplify(B*p^n + Bc*pc^n,'Steps', 200)



syms w real


s = exp(i*w);
sc = conj(s);

H0 = (s-z2)*(s-z1)*(sc-z1)*(sc-z2);


H = simplify(H0,'Steps',100)

Hd0 = (s-p)*(s-pc)*(sc-pc)*(sc-p);

Hd = simplify(Hd0,'Steps',150)
