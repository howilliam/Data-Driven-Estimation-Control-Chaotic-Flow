% POD coefficients plot
for j = 1:5
    switch j
        case 1 
            i = 1:10;
        case 2 
            i = 11:20;
        case 3
            i = 21:30;
        case 4
            i = 31:40;
        case 5
            i = 41:50;
    end   
    
    figure('name',['POD coefficients (modes ',num2str(i(1)),'-',num2str(i(end)),')'])
    for Nplot = 1:10
        subplot(5,2,Nplot);
        % Training dataset
        p1 = plot(Time(1:K_training),a(i(Nplot),1:K_training),'b-','Linewidth',1.5); hold on;
        p2 = plot(Time(1:K_training),a_e(i(Nplot),1:K_training),'r-','Linewidth',1.5); hold on;
        
        % Validation dataset
        p3 = plot(Time(K_training+1:end),a(i(Nplot),K_training+1:end),'b-','Linewidth',1.5); hold on;
        p4 = plot(Time(K_training+1:end),a_e(i(Nplot),K_training+1:end),'r-','Linewidth',1.5); hold on;
        
        plot(Time(K_training)*ones(1,2),[1.1*min(a(i(Nplot),:)), 1.1*max(a(i(Nplot),:))],'k-','Linewidth',2); % vertical line separating datasets
        if i == 1
            legend([p1,p2],'True','Estimate','Interpreter','Latex','FontSize',14,'Orientation','Horizontal');
        end
        ylabel(['${a}_{' num2str(m(i(Nplot))) '}(t)$'],'Interpreter','Latex','FontSize',18);
        xlim([Time(1) Time(end)]); ylim([min(a(i(Nplot),:)) max(a(i(Nplot),:))]);
        title(['FIT [\%]: ',num2str(round(FIT_a_training(i(Nplot)),2)),'\% (Training), ',...
            num2str(round(FIT_a_validation(i(Nplot)),2)),'\% (Validation)',...
            ],'Interpreter','Latex','FontSize',16);
        
        if Nplot == 9 || Nplot == 10
            xlabel('Time','interpreter','latex','FontSize',14);
        end
        set(gca,'linewidth',1);
        set(gca,'Fontsize',14);
        set(gca,'Ticklength',[0 0]);
    end
end