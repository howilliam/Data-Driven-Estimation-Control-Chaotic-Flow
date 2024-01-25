%%
% if type ==1
%     fileextension=['_' num2str(length(sensor_point)) '_sensors'];
% elseif type ==2
%     fileextension= ['optimal_' num2str(r) '_sensors'];
% end
% POD coefficients plot
for j = 1:2
    
    switch j
        case 1 
            i = 1:10; %%%%%%%%%%%%%%PILAS ACA 10
        case 2 
            i = 11:20;
        case 3
            i = 21:30;
        case 4
            i = 31:40;
        case 5
            i = 41:50;
    end   
    
    figure('name',['POD coefficients (modes ',num2str(i(1)),'-',num2str(i(end)),')'],'Units','normalized','Position',[0 0 1 1])
    
%     if i(1) == 11
%         var=10;
%     else 
%         var=10; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PILAS ACA 10
%     end
    
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
%     saveas(gcf,[pathvid 'Modes_' num2str(j) fileextension],'epsc')
%     saveas(gcf,[pathvid 'Modes' num2str(j) fileextension])
end

%% export plots

%%

% 
% % FIT [%] modes
% FIT_a_SID = zeros(length(m),1); FIT_a_KALMAN = FIT_a_SID; % pre-allocation
% for i = 1:length(m)
%     % SID
%     FIT_a_SID(i) = 100*(1-((norm(a(i,K_training+1:end) - a_SID(i,1+K_training:end)))/(norm(a(i,1+K_training:end) - mean(a(i,1+K_training:end))))));
%     
%     % Validation dataset KALMAN
%     FIT_a_KALMAN(i) = 100*(1-((norm(a(i,1+K_training:end) - a_e(i,1+K_training:end)))/(norm(a(i,1+K_training:end) - mean(a(i,K_training+1:end))))));
% end
% FIT_a = [FIT_a_SID, FIT_a_KALMAN];  
% 
% 
% 
% 
% if type ==1
%     fileextension=['_' num2str(length(sensor_point)) '_sensors'];
% elseif type ==2
%     fileextension= ['optimal_' num2str(r) '_sensors'];
% end
% % POD coefficients plot
% for j = 1:1
%     
%     switch j
%         case 1 
%             i = 1:10;
%         case 2 
%             i = 11:20;
%         case 3
%             i = 21:25;
% %         case 4
% %             i = 31:40;
% %         case 5
% %             i = 41:50;
%     end   
%     
%     figure('name',['POD coefficients (modes ',num2str(i(1)),'-',num2str(i(end)),')'],'Units','normalized','Position',[0 0 1 1])
%     
%     if i(1) == 11
%         var=10;
%     else 
%         var=10; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
%     
%     for Nplot = 1:var
%         subplot(5,2,Nplot);
%         % Training dataset
%         p1 = plot(Time(1:K_training),a(i(Nplot),1:K_training),'b-','Linewidth',1.5); hold on;
%         p2 = plot(Time(1:K_training),a_e(i(Nplot),1:K_training),'r--','Linewidth',1.5); hold on;
%         p5 = plot(Time(1:K_training),a_SID(i(Nplot),1:K_training),'k:','Linewidth',1.5); hold on;
%         % Validation dataset
%         p3 = plot(Time(K_training+1:end),a(i(Nplot),K_training+1:end),'b-','Linewidth',1.5); hold on;
%         p4 = plot(Time(K_training+1:end),a_e(i(Nplot),K_training+1:end),'r--','Linewidth',1.5); hold on;
%         p6 = plot(Time(K_training+1:end),a_SID(i(Nplot),K_training+1:end),'k:','Linewidth',1.5); hold on;
%         
%         plot(Time(K_training)*ones(1,2),[1.1*min(a(i(Nplot),:)), 1.1*max(a(i(Nplot),:))],'k-','Linewidth',2); % vertical line separating datasets
%         if Nplot == 1
%             legend([p3,p4,p6],'True','Kalman estimate','SID Estimate','Interpreter','Latex','FontSize',14,'Orientation','Horizontal');
%         end
%         ylabel(['${a}_{' num2str(m(i(Nplot))) '}(t)$'],'Interpreter','Latex','FontSize',18);
%         xlim([Time(1) Time(end)]); ylim([min(a_SID(i(Nplot),:)) max(a_SID(i(Nplot),:))]); %%"(££)(&")&)"(&)£"&)"(&")(&"
%         title(['FIT [\%]: ',num2str(round(FIT_a_SID(i(Nplot)),2)),'\% (SID), ',...
%             num2str(round(FIT_a_KALMAN(i(Nplot)),2)),'\% (Kalman)',...
%             ],'Interpreter','Latex','FontSize',16);
%         
%         if Nplot == var-1 || Nplot == var
%             xlabel('Time','interpreter','latex','FontSize',14);
%         end
%         set(gca,'linewidth',1);
%         set(gca,'Fontsize',14);
%         set(gca,'Ticklength',[0 0]);
%     end
% %     saveas(gcf,[pathvid 'Modes_' num2str(j) fileextension],'epsc')
% %     saveas(gcf,[pathvid 'Modes' num2str(j) fileextension])
% end
