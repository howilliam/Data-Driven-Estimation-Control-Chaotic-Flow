%% POD MODES SCRIPT 
% This is the second of a series of codes 

%% LOAD DATASETS 
%path='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\HPC data\DNS DATA (mat files)\';
path='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\';
pathPOD='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\';
load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\cyl_geo.mat')
load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\fishpoints_ORIG.mat')
load([path 'Griddata_HPC_200500_002TR20']);

%% DATASET size

x_length=[2 14];
y_length=[7 13]; % Domain coordinates for sensor grid in star-ccm
dx= 0.1; dy=0.1;
dt=0.01; %temporal resolution of CFD simulation
r=40; % rate at which the data is extracted
T=100020*dt:dt*r:250220*dt;% flow time
n=length(x_plotting);
K=length(T); %number of snapshots

%% Mean flow
Y_mean = mean([U;V],2);
% CAN INCLUDE A CHECK IF ALREADY CALCULATED

%% Snapshot
%crate snaphot matrix of velocity fluctuations

Y=[U;V] - Y_mean; 


% POD Calculation

    % Singular value decomposition (economy mode)
    [PHI,SIGMA]=svd(sqrt(dx*dy)*Y, 'econ'); 
    %phi is the matrx of singular vectors & sigma is the matrix of singular
    %values
    
    %Eigenvalues (energy of POD Modes)
    lambda= diag(SIGMA).^2*K;
    
    %Eigenvectors (spatial POD modes)
    phi = (1/sqrt(dx*dy) )*PHI;
    
    % Divergence of spatial POD modes (perform check for the first POD mode)
    div_check = sum(sum(divergence(xq,yq,griddata(x_plotting,y_plotting,phi(1:n,1),xq,yq),...
        griddata(x_plotting,y_plotting,phi(n+1:end,1),xq,yq))));
    
    % Orthonormality of eigenvectors
    orthonormal_check = zeros(10,10); % pre-allocation    
    for i = 1:10
        for j = 1:10
            orthonormal_check(i,j) = sum((phi(1:n,i).*phi(1:n,j) + ...
                phi(n+1:2*n,i).*phi(n+1:2*n,j))*dx*dy);
        end
    end
    
    %POD temporal coefficients 
    a= (Y'*sqrt(dx*dy)*PHI)';
    
  
        % Orthogonality of POD coefficients
    orthogonality_check = zeros(10,10); % pre-allocation
    for i = 1:10
        for j = 1:10
            orthogonality_check(i,j) = dot(a(i,:),a(j,:));
        end
    end

%% Energy of POD Modes 
MODES=5;
limit=25;
eigenvals_POD = figure('Name','Energy spectrum of POD modes'); 
subplot(1,2,1);
plot(1:MODES,lambda(1:MODES)/sum(lambda),'r-o','Linewidth',1,'Markersize',4,'MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
plot(MODES+1:length(lambda),lambda(MODES+1:end)/sum(lambda),'-o','Linewidth',1,'Color',[0.5 0.5 0.5],'Markersize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
xlabel('Mode $i$','fontsize',14,'interpreter','latex'); 
ylabel('$\lambda_{i}/\sum_{i=1}^{K}\lambda_i$','fontsize',14,'interpreter','latex'); 
grid on; xlim([1 limit]); xticks([0:5:limit]); ylim([0 0.5]);
yBox = [0, 0, 0.5, 0.5, 0];
xBox = [1, MODES, MODES, 1, 1];
retained = patch(xBox, yBox, 'red', 'FaceColor', 'red', 'EdgeColor','none','FaceAlpha', 0.2); 



yBox = [0, 0, 100, 100, 0];
xBox = [5, 10, 10, 5, 5];

retained2= patch(xBox, yBox, 'blue', 'FaceColor', 'blue', 'EdgeColor','none','FaceAlpha', 0.2); 
xBox = [10, 150, 150, 10, 10];
patch(xBox, yBox, 'black', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none','FaceAlpha', 0.2);
set(gca,'Fontsize',14,'TickLabelInterpreter','latex'); 


subplot(1,2,2);
plot(1:MODES,100*cumsum(lambda(1:MODES)/sum(lambda)),'r-o','Linewidth',1,'Markersize',4,'MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
energy_truncated = cumsum(lambda(1:MODES)/sum(lambda));
plot(MODES+1:length(lambda),100*(energy_truncated(end)+cumsum(lambda(MODES+1:end)/sum(lambda))),'-o','Linewidth',1,'Color',[0.5 0.5 0.5],'Markersize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
xlabel('Mode $i$','fontsize',14,'interpreter','latex'); 
ylabel('Cumulative Energy fraction (\%)','fontsize',14,'interpreter','latex'); 
grid on; xlim([1 limit]); xticks([0:5:limit]); ylim([0 100]);
yBox = [0, 0, 100, 100, 0];
xBox = [1, MODES, MODES, 1, 1];
patch(xBox, yBox, 'red', 'FaceColor', 'red', 'EdgeColor','none','FaceAlpha', 0.2); 
% xBox = [MODES, 150, 150, MODES, MODES];
% patch(xBox, yBox, 'black', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none','FaceAlpha', 0.2);


yBox = [0, 0, 100, 100, 0];
xBox = [5, 10, 10, 5, 5];
patch(xBox, yBox, 'blue', 'FaceColor', 'blue', 'EdgeColor','none','FaceAlpha', 0.2); 
xBox = [10, 150, 150, 10, 10];
truncated=patch(xBox, yBox, 'black', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none','FaceAlpha', 0.2);

legend([retained, retained2 ,truncated],'Reduced-order model','Retained for tranformation','Truncated modes','Fontsize',14,'interpreter','latex'); 
set(gca,'Fontsize',14,'TickLabelInterpreter','latex'); 
 %% POD spatial modes
PHIu = figure('Name','POD spatial modes - u velocity');
for i = 1:16
    subplot(4,4,i);
    zi = griddata(x_plotting,y_plotting,phi(1:n,i),xq,yq); 
    contourf(xq,yq,zi,MODES,'Edgecolor','flat'); hold on;  
    shading interp; axis equal; view(2); colormap(othercolor('RdBu10')); c = colorbar('eastoutside');
    set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')));
    xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
    xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
    title(['$\phi_{',num2str(i),'}$'],'fontsize',16,'interpreter','latex'); 
    set(gca,'linewidth',1);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontsize',14); 
    fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
end 

PHIv = figure('Name','POD spatial modes - v velocity');
for i = 1:16
    subplot(4,4,i); 
    zi = griddata(x_plotting,y_plotting,phi(n+1:end,i),xq,yq); 
    contourf(xq,yq,zi,MODES,'Edgecolor','flat'); hold on;
    shading interp; axis equal; view(2); colormap(othercolor('RdBu10')); c = colorbar('eastoutside');
    xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]); 
    xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
    set(c, 'YTickLabel', cellstr(num2str(reshape(get(c, 'YTick'),[],1),'%0.3f')));
    title(['$\phi_{',num2str(i),'}$'],'fontsize',16,'interpreter','latex'); 
    set(gca,'linewidth',1);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontsize',14); 
    fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;
    fill(x_cyl,y_cyl,'k');
end 

%% POD temporal coefficients
title_fig = ["POD temporal coefficients 1-5","POD temporal coefficients 6-10","POD temporal coefficients 11-15",...
    "POD temporal coefficients 16-20","POD temporal coefficients 21-25","POD temporal coefficients 26-30",...
    "POD temporal coefficients 31-35","POD temporal coefficients 36-40","POD temporal coefficients 41-45",...
    "POD temporal coefficients 46-50"];
for j = 1:2
    switch j
        case 1
            i = 1:5;
        case 2
            i = 6:10;
%         case 3
%             i = 11:15;
%         case 4
%             i = 16:20;
%         case 5
%             i = 21:25;
%         case 6
%             i = 26:30;
%         case 7
%             i = 31:35;
%         case 8
%             i = 36:40;
%         case 9
%             i = 41:45;
%         case 10
%             i = 46:50;
    end
    figure('Name',title_fig(j));
    
    for n = 1:5
        % Time domain
        subplot(5,2,2*n-1);
        plot(T,a(i(n),:),'k-','Linewidth',1);
        grid on; xlim([T(1) T(end)]);
        ylabel(['$a_{',num2str(i(n)),'}(t)$'],'fontsize',14,'interpreter','latex');
        if n == 5
            xlabel('Time','fontsize',14,'interpreter','latex');
        end
        set(gca,'Fontsize',12);
        
        % Frequency domain
        subplot(5,2,2*n);
        %___________________STROUHAL NUMBER = f*U/D; 
        u=1;
        D=2;
        [PSD] = spectrum_analyser([T',a(i(n),:)'],0)*D/u;
        s = plot(PSD(1,:),PSD(2,:),'k-','Linewidth',1); hold on;
        grid on; xlim([0 0.5]);
        ylabel(['$|a_{',num2str(i(n)),'}(St)|$'],'fontsize',14,'interpreter','latex');
        
        % insert datatip
        [~,indx] = max(PSD(2,:));
        dt = datatip(s,'DataIndex',indx,'interpreter','latex','fontsize',10);
        s.DataTipTemplate.DataTipRows(1).Label = '$St$'; s.DataTipTemplate.DataTipRows(1).Format = '%0.2f';
        s.DataTipTemplate.DataTipRows(2).Label = '$|a(St)|$'; s.DataTipTemplate.DataTipRows(2).Format = '%0.2f';
        if n == 5
            xlabel('$St$','fontsize',14,'interpreter','latex');
        end
        set(gca,'Fontsize',12);
    end
end

%% Fluctuations reconstruction from ROM
n = length(x_plotting); % number of grid points (excluding points inside the cylinders)
DNSvsROM = figure('Name','DNS vs. ROM velocity magnitude');
snapshot = find(round(T,4) == 200020*0.002); % choose snapshot to compare
Magnitude_DNS = sqrt((Y(1:n,snapshot) + Y_mean(1:n)).^2 + (Y(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % DNS

modes = [1,2,5,10]; % investigate different ROM sizes
for j = 1:length(modes)
    
    m = 1:modes(j); % modes retained by ROM
    Y_ROM = zeros(size(Y_mean)); % pre-allocation
    for i = 1:length(m)
        Y_ROM = Y_ROM + phi(:,m(i))*a(m(i),:); 
    end
    Magnitude_ROM = sqrt((Y_ROM(1:n,snapshot) + Y_mean(1:n)).^2 + (Y_ROM(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % ROM
    
    DNS_contours = subplot(4,3,[1:3:10]);
    zi = griddata(x_plotting,y_plotting,Magnitude_DNS,xq,yq);
    contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
    shading interp; axis equal; view(2); 
    colormap(DNS_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
    xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
    xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
    title('DNS','fontsize',14,'interpreter','latex');
    set(gca,'linewidth',1);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontsize',12); 
    fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
    
    ROM_contours = subplot(4,3,3*j-1); 
    zi = griddata(x_plotting,y_plotting,Magnitude_ROM,xq,yq);
    contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
    shading interp; axis equal; view(2); 
    colormap(ROM_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
    xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
    xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
    title(['ROM (',num2str(modes(j)),' modes)'],'fontsize',14,'interpreter','latex');
    set(gca,'linewidth',1);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontsize',12); 
    fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
    
    FITTING = subplot(4,3,3*j);
    FIT = 100*(1 - abs(Magnitude_ROM - Magnitude_DNS)./abs(Magnitude_DNS));
    zi = griddata(x_plotting,y_plotting,FIT,xq,yq);
    zi = max(zi,0); % remove all negative fit values
    contourf(xq,yq,zi,80,'Edgecolor','flat'); hold on;
    shading interp; axis equal; view(2); colormap(FITTING,jet); c = colorbar('eastoutside'); 
    caxis([0 100]); colorbar('Ticks',[0:25:100]);
    xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
    xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
    title('FIT [\%]','fontsize',14,'interpreter','latex');
    set(gca,'linewidth',1);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontsize',12); 
    fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
end


%% Save variables
% % Mean flow 
% save('Y_mean.mat','Y_mean'); 
% 
% % Snapshot matrix of fluctuations
% save('Y.mat','Y');
% 
% % Eigenvalues 
% save('lambda.mat','lambda');
% 
% % Eigenvectors (only ROM modes)
%phi = phi(:,1:50);
% save('phi.mat','phi'); 
% 
% % Temporal coefficients
% save('a.mat','a');

%save([pathPOD 'POD_HPC_100500_002TR20'],'a','phi','Y_mean','Y','lambda')