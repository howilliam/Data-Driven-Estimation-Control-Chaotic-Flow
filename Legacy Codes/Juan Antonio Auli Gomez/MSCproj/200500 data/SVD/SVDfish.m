%% SVD fish sensors
function [a_P,Y_P,Y_P_mean,lam_p,Phi_P]=SVDfish(P_fish,snapshot,MODES,X_fish,Y_fish,TYPE)
%%Temporal definitions
dt=0.002; %temporal resolution of data extraction
r=40; % rate at which the data is extracted
T=100020*dt:dt*r:250220*dt;% flow time
K=length(T); %number of snapshots
path='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\100s data\';
%% Data
    
    % Create snapshots
    Y_P_mean= mean(P_fish,2);
    Y_P= P_fish - Y_P_mean;
    
    %Perform Singule value decomposition
    [Phi_P,Sig_P]=svd(Y_P,'econ');
        
        % Calculate the Eigenvalues
        lam_p= diag(Sig_P) .^2*K;
        
        % eigenvectors 
        
        %Temporal coefficient 
        a_P=(Y_P'*Phi_P)'; %NB (Y_P'*Phi_P)' = Phi_P'*Y_P
        
        %Orthogonality Check 
            orthogonality_check = zeros(10,10); % pre-allocation
            for i = 1:10
                for j = 1:10
                    orthogonality_check(i,j) = dot(a_P(i,:),a_P(j,:));
                end
            end

        



%% 'ENERGY' OF MODES (SVD -)
limit=25;
%load([path 'fish_sorted_TR20_30June.mat'], 'X_fish', 'Y_fish');
eigenvals_SVDfish = figure('Name','Energy spectrum of POD modes'); 
subplot(1,2,1);
plot(1:MODES,lam_p(1:MODES)/sum(lam_p),'r-o','Linewidth',1,'Markersize',4,'MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
plot(MODES+1:length(lam_p),lam_p(MODES+1:end)/sum(lam_p),'-o','Linewidth',1,'Color',[0.5 0.5 0.5],'Markersize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
xlabel('Mode $i$','fontsize',14,'interpreter','latex'); 
ylabel('$\lambda_{i}/\sum_{i=1}^{K}\lambda_i$','fontsize',14,'interpreter','latex'); 
grid on; xlim([1 limit]); xticks([0:5:limit]); ylim([0 0.5]);
yBox = [0, 0, 0.5, 0.5, 0];
xBox = [1, MODES, MODES, 1, 1];
retained = patch(xBox, yBox, 'red', 'FaceColor', 'red', 'EdgeColor','none','FaceAlpha', 0.2); 
xBox = [MODES+5, 30, 30, MODES+5, MODES+5];
truncated = patch(xBox, yBox, 'black', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none','FaceAlpha', 0.2);


yBox = [0, 0, 100, 100, 0];
xBox = [5, 10, 10, 5, 5];
retained2= patch(xBox, yBox, 'blue', 'FaceColor', 'blue', 'EdgeColor','none','FaceAlpha', 0.2); 

legend([retained, retained2,truncated],'Reduced-order model','Retained for tranformation','Truncated modes','Fontsize',14,'interpreter','latex'); 
set(gca,'Fontsize',14,'TickLabelInterpreter','latex'); 

subplot(1,2,2);
plot(1:MODES,100*cumsum(lam_p(1:MODES)/sum(lam_p)),'r-o','Linewidth',1,'Markersize',4,'MarkerFaceColor','r','MarkerEdgeColor','r'); hold on;
energy_truncated = cumsum(lam_p(1:MODES)/sum(lam_p));
plot(MODES+1:length(lam_p),100*(energy_truncated(end)+cumsum(lam_p(MODES+1:end)/sum(lam_p))),'-o','Linewidth',1,'Color',[0.5 0.5 0.5],'Markersize',4,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
xlabel('Mode $i$','fontsize',14,'interpreter','latex'); 
ylabel('Cumulative Energy fraction (\%)','fontsize',14,'interpreter','latex'); 
grid on; xlim([1 limit]); xticks([0:5:limit]); ylim([0 100]);
yBox = [0, 0, 100, 100, 0];
xBox = [1, MODES, MODES, 1, 1];
patch(xBox, yBox, 'red', 'FaceColor', 'red', 'EdgeColor','none','FaceAlpha', 0.2); 
xBox = [MODES+5, 30, 30, MODES+5, MODES+5];
patch(xBox, yBox, 'black', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none','FaceAlpha', 0.2);
yBox = [0, 0, 100, 100, 0];
xBox = [5, 10, 10, 5, 5];
retained2= patch(xBox, yBox, 'blue', 'FaceColor', 'blue', 'EdgeColor','none','FaceAlpha', 0.2); 


set(gca,'Fontsize',14,'TickLabelInterpreter','latex'); 
%% Reduced order model 

p=length(X_fish); % number of sensors in fish;

    %Comparison of reconsruction
    
    pointsize = 20;  %<====================================================
        FISH_MODES = [3,5,10,15]; % investigate different ROM sizes
        var=figure('name','DNS vs ROM and FIT %');
        FISH_DNS_Pressure = subplot(4,3,[7 10]);
        scatter(X_fish, Y_fish, pointsize, Y_P(:,snapshot)+Y_P_mean,'o','filled'); axis equal
        colorbar('eastoutside');cb=colormap(othercolor('RdBu10')); title(['DNS ' TYPE ' along fish'],'fontsize',14,'interpreter','latex')
        %axis off; 
        
%         DNS_contours = subplot(4,3,[1 4]);
%         
%         load([path 'Griddata_TR20_june24.mat'],'xq','yq','x_plotting','y_plotting');
%         load([path 'POD_TR_20_june24.mat'],'Y','Y_mean');
%         n = length(x_plotting);
%         Magnitude_DNS = sqrt((Y(1:n,snapshot) + Y_mean(1:n)).^2 + (Y(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % DNS
%         zi = griddata(x_plotting,y_plotting,Magnitude_DNS,xq,yq);
%         contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
%         shading interp; axis equal; view(2); 
%         colormap(DNS_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([min(Magnitude_DNS) max(Magnitude_DNS)]);
%         x_length=[2 14];
%         y_length=[7 13];
%         xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
%         xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
%         title('DNS velocity magnitue','fontsize',14,'interpreter','latex');
%         set(gca,'linewidth',1);
%         set(gca,'TickLength',[0 0]);
%         set(gca,'Fontsize',12); load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\cyl_geo.mat');
%         %load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\fishpoints_ORIG.mat');
%         load('/home/jaa21/Downloads/Project/additionalfiles/fishpoints_ORIG.mat')
%         fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');;hold on;fill(x_cyl,y_cyl,'k');
        
    
    for j = 1:length(FISH_MODES)
        %Reduced order model (in the different modes to be compared)
        FISH_ROM_PRESSURE = subplot(4,3,3*j-1);
        Y_P_ROM= Phi_P(:,1:FISH_MODES(j))*a_P(1:FISH_MODES(j),:);
        scatter(X_fish, Y_fish, pointsize, Y_P_ROM(:,snapshot)+Y_P_mean,'o','filled');axis equal
        colormap(FISH_ROM_PRESSURE,othercolor('RdBu10')); colorbar('eastoutside'); caxis([min(Y_P(:,snapshot)+Y_P_mean) max(Y_P(:,snapshot)+Y_P_mean)]);
        title(['ROM (' num2str(FISH_MODES(j)) ' modes)'],'fontsize',14,'interpreter','latex');%axis off


        % FIT 
        FIT=100* (1-abs((Y_P(:,snapshot) - Y_P_ROM(:,snapshot)))./abs((Y_P(:,snapshot) +Y_P_mean)) );
        FITTING = subplot(4,3,3*j);
        scatter(X_fish, Y_fish, pointsize, FIT(),'o','filled'); axis equal
        colorbar('eastoutside'); caxis([0 100]); colormap(FITTING,jet);title('FIT [\%]','fontsize',14,'interpreter','latex');%axis off
    end
    
    %%


        
%% MAX POD (pressure)
% clear POD_peaks;
% [POD_peaks] = maxPOD_fish(Phi_P,9);
% max_POD_Pressure=figure('name','Max POD distribution along fish')
% plot([s_geo s_geo],[h_geo_L h_geo_u]); hold on;
% 
% xx=[];
% yy=[];
% for i=1:length(POD_peaks)
% %scatter(X_fish(POD_peaks(i)),Y_fish(POD_peaks(i)),'ro');
% xx=[xx;X_fish(POD_peaks(i))];
% yy=[yy;Y_fish(POD_peaks(i))];
% 
% hold on; axis equal
% 
% end 
% scatter(xx,yy,'ro')
% labels=1:1:length(POD_peaks);
% labelpoints(xx,yy,labels)

