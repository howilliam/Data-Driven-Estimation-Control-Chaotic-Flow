function [fig]=sensorplot(indecies,title_plot)



% load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\fish_sorted_HPC_200500_002TR20.mat','Y_fish','X_fish');% contains data order as used in codes
% fig=figure;
% plot(X_fish([1:2:500 500]),Y_fish([1:2:500 500]),'k', 'linewidth', 1.3); hold on;
% plot(X_fish([1 2:2:500]),Y_fish([1 2:2:500]),'k', 'linewidth', 1.3); axis equal
% ylim([9.5 10.5]); 
% xlim([7.5 12.5])

title(title_plot,'interpreter','latex','fontsize',15)

s=scatter(X_fish(indecies),Y_fish(indecies),'ro');

% labelpoints(X_fish(indecies),Y_fish(indecies),[1:length(indecies)],'Color','k','buffer',-0.35,'FontWeight','bold','interpreter','latex','FontSize',14); %,'Color','r'
fish(indecies(i)),Y_fish(indecies(i)),[i],'Color','k','buffer',0.10,'position','C','FontWeight','bold','interpreter','latex','FontSize',14);
end
set(gca,'XTick',[], 'YTick', [])


end