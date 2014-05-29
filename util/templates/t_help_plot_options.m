% text(ix,iy,num2str( Sn(iy,ix)),'Color','w','Interpreter','latex','FontSize',14,'BackgroundColor',[0.2 0.2 0.2],'HorizontalAlignment','Center')
% plot(scale*z(1,:) + xind1,   scale*z(2,:) + yind1,'linewidth',2,'color','yellow')
% plot(xind1, yind1, 'wo',  'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4)   
% text(ix,iy,num2str(Sn(iy,ix)),'Color','w','Interpreter','latex','FontSize',14,'BackgroundColor',[0.2 0.2 0.2],'HorizontalAlignment','Center')
% plot(xind, yind, 'o', 'MarkerFaceColor','g')

% title('80%','fontsize',12,'fontweight','b')
% set(gca,'fontsize',12,'fontweight','b', 'XTick', (0:100:1000) )

%% using a single cell variable for the setting
% plot_option={'linewidth',4};
% plot(rcoord(:,1),plot_option{:}),hold on
% plot(rcoord(:,2),'r',plot_option{:})
% plot(rcoord(:,3),'g',plot_option{:})
%% setting axes options
% axesOption={'XTick',[],'YTick',[],'linewidth',3};
% figure,imexpl,spy(L_brute),set(gca,axesOption{:}),xlabel(''),box on
%% saving colorbar images
% CAXIS=0.02*[-1 1];
% load('elasticnet_admm_sahd.mat')
% imcovvl(v2mat,CAXIS)
% cbar=colorbar('location','west','fontsize',25','fontweight','b','YTick',(-.02:.01:.02))
% % export_fig(cbar, 'colorbar.pdf','-transparent',true,'-nocrop',true
% 
% % http://stackoverflow.com/questions/13648149/define-dimensions-of-colorbar-in-matlab
% x1=get(gca,'position');
% x=get(cbar,'Position');
% x(3)=0.01;
% set(cbar,'Position',x)
% set(gca,'position',x1)
% export_fig(cbar, 'colorbar.png','-transparent',true,'-r150');