clear; close all;
w0 = pi/4;
r = 0.8;
r2 = 0.5;
p = r*exp(j*w0);
p2 = r2*exp(j*w0);
w = 0:.01:pi;
z = exp(j*w);

H1 = 1-p*z;
H2 = z-p;
H = (1-p*z)./(z-p);



T = -w - atan( (r*sin(w+w0))./(1-r*cos(w+w0)) ) - atan( (r*sin(w-w0))./(1-r*cos(w-w0)) );
T2 = -w + atan( (r2^2*sin(2*w)-2*r2*sin(w)*cos(w0))./(1-2*r2*cos(w)*sin(w0)+r2^2*cos(2*w)) );
T3 = -w - 2*atan( (r*sin(w-w0))./(1-r*cos(w-w0)) );
figure;
plot(w/pi,angle(H));hold on;
%plot(w/pi,T3,'LineWidth',2); hold on; grid on;
plot(w/pi,T,'LineWidth',2);
%plot(w/pi,unwrap(T2),'LineWidth',2);
plot([w0/pi,w0/pi], [-3.5,0],'k--')

xlabel('\omega/\pi'); ylabel('\Phi(\omega)');
legend('r = 0.7', 'r = 0.9')
title('First Order Complex All-Pass Phase Response')

% 
% figure,plot(w,angle(H1));hold on;
% plot(w,unwrap(atan2(-r*sin(delta),(1-r*cos(delta)))))
% plot(w,angle(H2));plot(w,(atan2((sin(w)-r*sin(w0)),(cos(w)-r*sin(w0)))));