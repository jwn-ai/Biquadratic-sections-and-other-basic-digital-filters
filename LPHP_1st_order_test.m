clear; close all;
a = 0.5;
a2 = 0.9;
a3 = 0.1;
w = .0001:.01:pi;
z = exp(j*w);

Hlp0 = (1-a)/2*(1+z.^(-1))./(1-a*z.^(-1));
Hhp0 = (1+a)/2*(1-z.^(-1))./(1-a*z.^(-1));

Tlp0 = @(a) -atan2( (1+a)*sin(w),(1-a)*(1+cos(w)));
Thp0 = @(a) atan2( (1-a)*sin(w),(1+a)*(1-cos(w)));

figure,plot(w/pi,angle(Hlp0));
figure,plot(w/pi,abs(Hlp0));
Hlp = @(a,b) (1-a)/(1+b) * sqrt((1+2*b*cos(w)+b^2)./(1-2*a*cos(w)+a^2));
Hhp = @(a,b) (1+a)/(1+b) * sqrt((1-2*b*cos(w)+b^2)./(1-2*a*cos(w)+a^2));

H1 = Hlp(a,1);
H2 = Hlp(a2,1);
H3 = Hlp(a3,1);

G1 = Hhp(a,1);
G2 = Hhp(a2,1);
G3 = Hhp(a3,1);


figure;
plot(w/pi,H1,'k-','LineWidth',2);hold on;
plot(w/pi,H2,'k--','LineWidth',2);grid on;
plot(w/pi,H3,'k-.','LineWidth',2);grid on;


xlabel('\omega/\pi'); ylabel('G(\omega)');
legend('a = 0.5', 'a = 0.9', 'a = 0.2')
title('First Order IIR LPF Magnitude')

figure;
%plot(w/pi,unwrap(angle(Hlp0),.5));hold on;
plot(w/pi,Tlp0(a),'k-','LineWidth',2);hold on;grid on;
plot(w/pi,Tlp0(a2),'k--','LineWidth',2);
plot(w/pi,Tlp0(a3),'k-.','LineWidth',2);

xlabel('\omega/\pi'); ylabel('\Phi(\omega)');
legend('a = 0.5', 'a = 0.9', 'a = 0.2')
title('First Order IIR LPF Phase')

figure;
plot(w/pi,G1,'k-','LineWidth',2);hold on;
plot(w/pi,G2,'k--','LineWidth',2);grid on;
plot(w/pi,G3,'k-.','LineWidth',2);grid on;

xlabel('\omega/\pi'); ylabel('G(\omega)');
legend('a = 0.5', 'a = 0.9', 'a = 0.2')
title('First Order IIR HPF Magnitude')

figure;
%plot(w/pi,unwrap(angle(Hlp0),.5));hold on;
plot(w/pi,Thp0(a),'k-','LineWidth',2);hold on;grid on;
plot(w/pi,Thp0(a2),'k--','LineWidth',2);
plot(w/pi,Thp0(a3),'k-.','LineWidth',2);

xlabel('\omega/\pi'); ylabel('\Phi(\omega)');
legend('a = 0.5', 'a = 0.9', 'a = 0.2')
title('First Order IIR HPF Phase')

