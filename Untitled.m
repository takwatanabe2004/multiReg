clear
purge

k=1:500;
% k=1:100:50e3;

figure,imexpl
subplot(311),tplot(k,log10(1./k))
subplot(312),tplot(k,log10(1./k.^2))
subplot(313),tplot(k,log10(1./k)),hold on,plot(k,log10(1./k.^2),'r')
% plot(k,log10(80./k.^2),'r')
% figure,imexpl
% subplot(211)

c=0.95;
figure,imexpl
tplot(k,log10(1./k)),hold on,plot(k,log10(1./k.^2),'r')
plot(k,log10(c.^k),'g')
legend('1/k','1/k^2','linear')