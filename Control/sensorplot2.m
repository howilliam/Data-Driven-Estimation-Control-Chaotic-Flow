function [fig]=sensorplot2(indecies,title_plot,p_modes)



load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\fish_sorted_HPC_200500_002TR20.mat','Y_fish','X_fish');% contains data order as used in codes
fig=figure;
plot(X_fish([1:2:500 500]),Y_fish([1:2:500 500]),'k', 'linewidth', 1.3); hold on;
plot(X_fish([1 2:2:500]),Y_fish([1 2:2:500]),'k', 'linewidth', 1.3); axis equal
ylim([9.5 10.5]); 
xlim([7.5 12.5])

title(title_plot,'interpreter','latex','fontsize',15)
if p_modes <= 5
s=scatter(X_fish(indecies),Y_fish(indecies),'ro');

else
    s=scatter(X_fish(indecies(1:5)),Y_fish(indecies(1:5)),'ro');
    s=scatter(X_fish(indecies(6:end)),Y_fish(indecies(6:end)),'bo');
end
% labelpoints(X_fish(indecies),Y_fish(indecies),[1:length(indecies)],'Color','k','buffer',-0.35,'FontWeight','bold','interpreter','latex','FontSize',14); %,'Color','r'
for i=1:p_modes
    if Y_fish(indecies(i))>10
        alph=1.01;
        beta=1;
    else
        alph=0.99;
        beta=1;
    end
    
    if X_fish(indecies(i))<8.08
        beta=0.98;
        if Y_fish(indecies(i))>10
        alph=1.001;
        else
            alph=0.999;
        end
    end
    
labelpoints(X_fish(indecies(i))*beta,Y_fish(indecies(i))*alph,[i],'Color','k','buffer',0.10,'position','C','FontWeight','bold','interpreter','latex','FontSize',14);
end
set(gca,'XTick',[], 'YTick', [])


end