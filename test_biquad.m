clear all;close all;
w0 = pi/4;
r = .9;
p = r*exp(j*w0);
pc = conj(p);
z1 = 1;
z2 = -1;

n = (0:1:50);

C = (p-pc)*p;

B = (p-z1)*(p-z2)/C;
rn = r.^n;
A = rn.*exp(j*w0*n);
Ac = conj(A);
Bc = conj(B);

z12= z1*z2;
zp12 = z1+z2;


delta = @(x) x==0;

f0 = delta(n)*z1*z2/r^2 + B.*A + Bc.*Ac;

w = -pi:.001:pi;


[A1,f1,argH1,H1] = biquad(n,z1,z2,r,pi/4,w);
[A2,f2,argH2,H2] = biquad(n,z1,z2,.75,pi/4,w);
figure(2);
stem(n,f1,'LineWidth',2,'Color','k');hold on;grid on;
%stem(n,f2,'LineWidth',2);



p =3/4;
b0 = 1;
h = b0*(n+1).*p.^n;
%figure;stem(n,h)

figure;
plot(w/pi,A1,'k-','LineWidth',2);hold on;grid on;
plot(w/pi,A2,'k--','LineWidth',2);
legend('r = 0.9','r = 0.75')
xlabel('\omega/\pi')
ylabel('|H(\omega)|')


figure;
plot(w/pi,angle(H1));hold on;grid on;
plot(w/pi,angle(H2))
plot(w/pi,argH1,'--')

function [Hmag,f1,argH,H] =  biquad(n,z1,z2,r,w0,w)
z12= z1*z2;
rn = r.^n;
zp12 = z1+z2;
delta = @(x) x==0;

f1 = delta(n)*z1*z2/r^2 + rn.*(sin(w0*(n+1)) - zp12/r*sin(n*w0) + z12/r^2*sin((n-1)*w0))/sin(w0);



Hmag = (((z1^2-2*z1*cos(w)+1).*(z2^2-2*z2*cos(w)+1))./...
    ((r^2-2*r*cos(w-w0)+1).*(r^2-2*r*cos(w+w0)+1))).^(1/2);
Hm = max(Hmag);
Hmag = Hmag/Hm;

argHd = atan2((sin(w)-r*sin(w0)),((cos(w)-r*cos(w0))))+...
    atan2((sin(w)+r*sin(w0)),((cos(w)-r*cos(w0))));
argHn = atan2(sin(2*w)-(z1+z2)*sin(w),cos(2*w)-(z1+z2)*cos(w)+z2*z1);
argH = unwrap(argHn - argHd);


p = r*exp(j*w0);
pc = conj(p);
s = exp(j*w);

H = ((s-z1).*(s-z2))./((s-p).*(s-pc));

arg1 = angle(s-p);
arg2 = angle(s-pc);
arg3 = angle((s-z1).*(s-z2));
%argH = -unwrap(arg3-arg1-arg2);
% 
% figure;
% plot(w,arg1+arg2);hold on;
% plot(w,arg3)
% plot(w,atan2(sin(2*w),cos(2*w)-1));hold on;
% plot(w,argHd)
% plot(w,unwrap(arg3-arg1-arg2))
% plot(w,argHn)
% plot(w,argHn-argHd)
% plot(w,unwrap(argHn-argHd))
% plot(w,argH)
end