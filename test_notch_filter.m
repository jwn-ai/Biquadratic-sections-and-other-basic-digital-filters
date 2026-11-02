
clear; close all;
w0 = pi/4;
r = 0.7;
r2 = 0.9;
p = r*exp(j*w0);
z1 = exp(j*w0);
pc = conj(p);
z1c = conj(z1);
p2 = r2*exp(j*w0);
p2c = conj(p2);
w = -pi:.001:pi;
z = exp(j*w);

%first order
G = abs(1+p)/2*(z-z1)./(z-p);
figure,plot(w,abs(G));

abs(1+p)/abs(1+z1)*abs(1-z1)/abs(1-p)

b01 = sqrt((1+2*r*cos(w0)+r^2)/(2+2*cos(w0)));
b02 = b01^2;

H0 = (z-z1)./(z-p);
H2 = b02*abs(((z-z1).*(z-z1c))./((z-p).*(z-pc)));
H1 = b01*abs((z-z1)./(z-p));
G2 = @(p) ((1+abs(p))/2)^2*abs(((z-z1).*(z-z1c))./((z-p).*(z-conj(p))));
G1 = @(p) ((1+abs(p))/2)*abs((z-z1)./(z-p));

T0 = angle(H0);
T = atan( (cos(w)-r*cos(w0))./(sin(w)-r*sin(w0)) )...
    - atan( (cos(w)-1*cos(w0))./(sin(w)-1*sin(w0)) );
T2 = atan( (cos(w)-r2*cos(w0))./(sin(w)-r2*sin(w0)) )...
    - atan( (cos(w)-1*cos(w0))./(sin(w)-1*sin(w0)) );

P2 = @(r) atan2( (cos(w)-1*cos(w0)),(sin(w)-1*sin(w0)) )...
    + atan2( (cos(w)-1*cos(w0)),(sin(w)+1*sin(w0)) )...
   - atan2( (cos(w)-r*cos(w0)),(sin(w)-r*sin(w0)) )...
    - atan2( (cos(w)-r*cos(w0)),(sin(w)+r*sin(w0)) );

figure;
plot(w,unwrap(P2(r)),'LineWidth',2);hold on; grid on;
plot(w,unwrap(P2(r2)),'LineWidth',2);
xlabel('\omega');ylabel('\angle H(\omega)')
legend('r = 0.7','r = 0.9')
title('2^{nd} Order Phase Response')


figure,plot(w,(T),'LineWidth',2);hold on;grid on;
plot(w,T2,'LineWidth',2);
xlabel('\omega');ylabel('\angle H(\omega)')
legend('r = 0.7','r = 0.9')
title('1st Order Phase Response')

figure;
plot(w,abs(H1),'LineWidth',2);hold on;grid;
plot(w,abs(H2),'LineWidth',2);
legend('First Order','Second Order')
xlabel('\omega');ylabel('|H(\omega)|')
title('DC Null Filter')
% 
% bw1 = x(1).Position(1)-x(2).Position(1)
% bw2 = x(3).Position(1)-x(4).Position(1)
% 
% bw2/bw1
B = sqrt(2);
B1 = 1-B;
rr = 0.1:.01:.99;
BW2 = (B1*rr.^2 + 2*rr +B1)./(rr.^2+2*rr*B1 +1);
BW1 = (2*rr)./(rr.^2+1);

w1 = acos(BW1);
w2 = acos(BW2);

rt = sec(.105166) -tan(.105166)

figure
plot(rr,w1,'LineWidth',2); grid on;hold on;
plot(rr,w2,'LineWidth',2);
xlabel('r');ylabel('Half-angle BW [rad]')
legend('First Order','Second Order')
title('Pole Magnitude vs. Bandwidth')



figure;
plot(w,G2(p),'LineWidth',2);hold on;grid;
plot(w,G2(p2),'LineWidth',2);
legend('r = 0.7','r = 0.9')
xlabel('\omega');ylabel('|H(\omega)|')
title('Second Order Notch Filter')


figure;
plot(w,G1(p),'LineWidth',2);hold on;grid;
plot(w,G1(p2),'LineWidth',2);
legend('r = 0.7','r = 0.9')
xlabel('\omega');ylabel('|H(\omega)|')
title('1^{st} Order Notch Filter')