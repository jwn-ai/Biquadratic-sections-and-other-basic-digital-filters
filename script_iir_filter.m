%function output = IIRFilter(z,filter)
clear all;
t = 0:99;
N = length(t);
t = t';
a = 1;
b = 0;
s = a*t+b;
alpha = 0.75;
sigma = 1;
Nruns = 1000;
for i=1:Nruns
    i
    z = s + sigma*randn(N,1);
    xhat = z(1);
    output.xhat(1) = xhat;
    shat = s(1);
    output.shat(1) = shat;
    for n=2:N; %frame index
        zn = z(n);
        xhat = zn*(1-alpha) + xhat*alpha;
        shat = s(n)*(1-alpha) + shat*alpha;
        output.xhat(n,1) = xhat;
        output.shat(n,1) = shat;
    end
    e = output.xhat - output.shat;
    e2 = output.xhat - s;
    e_avg(i,:) = e';
    e2_avg(i,:) = e2';
end
E = (1-alpha)./(1+alpha)*sigma^2 + alpha^2/(1-alpha)^2;
sigmax2 = (1-alpha)/(1+alpha)*sigma^2
e_avg = sum(e_avg.^2)/Nruns;
e2_avg = sum(e2_avg.^2)/Nruns;
figure,plot(sqrt(sigmax2)*ones(100,1),'r');
hold on;
plot(sqrt(E)*ones(100,1),'g');plot(sqrt(e_avg));
plot(sqrt(e2_avg),'m')
grid;
xlabel('samples');ylabel('error');
legend('/sigma_n^{2}','E','rms signal error','rms estimation error')


a = 0.01:.01:.5;
figure,plot(a,a.^2./(1-a).^2,'r');hold on; plot(a, (1-a)./(1+a) ,'b');
xlabel('\alpha');legend('bias','variance');