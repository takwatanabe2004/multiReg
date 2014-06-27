clear
purge
load Results_tak
Whos
figure,imexpb
subplot(121),imcov(ConnMean)
subplot(122),imcov(ConnVar)
% coord
% coord_used
figure,imexpb
subplot(121),imedge(edgemask),axis on
subplot(122),tplot(nodemask)
figure,imexpb
subplot(121),tplot(edgemask(163,:)'-nodemask)
subplot(122),tplot(edgemask(163,:))
% nodemask
isequal(edgemask(163,:)',nodemask)