function [fig]=sensorplot2(index_u,index_v,title_plot,p_modes)
% function [fig]=sensorplot2(index,title_plot,p_modes)


%load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\fish_sorted_HPC_200500_002TR20.mat','y_plotting','x_plotting');% contains data order as used in codes
load('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\Save Test\Griddata_HPC_200500_002TR20.mat','y_plotting','x_plotting');% contains data order as used in codes
fig=figure;
% plot(x_plotting([1:2:500 500]),y_plotting([1:2:500 500]),'k', 'linewidth', 1.3); hold on;
% plot(x_plotting([1 2:2:500]),y_plotting([1 2:2:500]),'k', 'linewidth', 1.3); axis equal
% ylim([9.5 10.5]); 
% xlim([7.5 12.5])

title(title_plot,'interpreter','latex','fontsize',15)
if p_modes <= 5
s=scatter(x_plotting(index_u),y_plotting(index_u),'ro');

else
    %s=scatter(x_plotting(index(1:5)),y_plotting(index(1:5)),'ro');
    s=scatter(x_plotting(index_u),y_plotting(index_u),'filled','rs');  %% u velocity index being plotted for sensors
    %s=scatter(x_plotting(index(6:end)),y_plotting(index(6:end)),'bo');
    hold on;
    index_v = index_v - length(x_plotting);
    s=scatter(x_plotting(index_v),y_plotting(index_v),'filled','bo');    %% v velocity index being plotted for sensors
    hold on;
    pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
    plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
    pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
    plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
    grid off; axis equal; view(2);
    xlim([-3 11]); ylim([-4 4]);
    xticks(0:2:round(14)); yticks(-4:2:4);
    
    set(gca,'linewidth',1);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontsize',14); 
    
end
% labelpoints(x_plotting(index),y_plotting(index),[1:length(index)],'Color','k','buffer',-0.35,'FontWeight','bold','interpreter','latex','FontSize',14); %,'Color','r'
% for i=1:p_modes
%     if y_plotting(index(i))>10
%         alph=1.01;
%         beta=1;
%     else
%         alph=0.99;
%         beta=1;
%     end
%     
%     if x_plotting(index(i))<8.08
%         beta=0.98;
%         if y_plotting(index(i))>10
%         alph=1.001;
%         else
%             alph=0.999;
%         end
%     end
%     
% labelpoints(x_plotting(index(i))*beta,y_plotting(index(i))*alph,[i],'Color','k','buffer',0.10,'position','C','FontWeight','bold','interpreter','latex','FontSize',14);
% end
% set(gca,'XTick',[], 'YTick', [])


end