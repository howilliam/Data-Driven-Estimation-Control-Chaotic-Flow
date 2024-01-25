
subplot(1,2,1);
        plot(TCD(1000:end),CD(1000:end),'k-','Linewidth',1);
        grid on; xlim([TCD(1000) TCD(end)]);
        ylabel(['$CD$'],'fontsize',14,'interpreter','latex');
        if n == 1
            xlabel('Time','fontsize',14,'interpreter','latex');
        end
        set(gca,'Fontsize',12);
        
        % Frequency domain
subplot(1,2,2);
        %___________________STROUHAL NUMBER = f*U/D; 
        u=1;
        D=2;
        [PSD_CD] = spectrum_analyser([TCD(2500:end),CD(2500:end)],0)*D/u;
        B = plot(PSD_CD(1,:),PSD_CD(2,:),'k-','Linewidth',1); hold on;
        grid on; xlim([0 0.5]);
        ylabel(['$|a_{',num2str(i(n)),'}(St)|$'],'fontsize',14,'interpreter','latex');
        
        % insert datatip
        [~,indx] = max(PSD_CD(2,:));
        dt = datatip(B,'DataIndex',indx,'interpreter','latex','fontsize',10);
        B.DataTipTemplate.DataTipRows(1).Label = '$St$'; B.DataTipTemplate.DataTipRows(1).Format = '%0.2f';
        B.DataTipTemplate.DataTipRows(2).Label = '$|a(St)|$'; B.DataTipTemplate.DataTipRows(2).Format = '%0.2f';
        if n == 5
            xlabel('$St$','fontsize',14,'interpreter','latex');
        end
        set(gca,'Fontsize',12)