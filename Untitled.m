clear
purge

t=linspace(-5,5,401);
lwid=2;
tau=0.5;
y = 2;
plot(t,tak_soft(t,tau),'linewidth',lwid), hold on,grid on
plot(t,y+tak_soft(t-y,tau),'r','linewidth',lwid)
y=-2;
plot(t,y+tak_soft(t-y,tau),'g','linewidth',lwid)
legend('soft','y=2','y=-2')
% tplott(t, abs(t-y))
% 
% tplott(t, abs(t+y))