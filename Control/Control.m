clc; close all;

clear; 
path   ='/home/wh219/FYP/Controller/Controller Save/'; % POD location
path2 = '/home/wh219/FYP/Controller/Controller Save/'; %grid data loction 
pathSID='/home/wh219/FYP/Controller/Controller Save/'; %location of the saved system identification
%% pathKalman='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\Kalman\';
load([path 'POD_Control_v3_9000.mat'],'a','phi','Y_mean','Y');
load([path2 'Griddata_Control_v3_9000.mat'],'xq','yq','x_plotting','y_plotting')
    
load('/home/wh219/FYP/Codes/MSCproj/200500 data/Control/control_stuff.mat');
load('/home/wh219/FYP/Codes/MSCproj/200500 data/Control/u_signal2.mat');

%% TEMPORAL DEFINITIONS
dt=0.01; %temporal resolution of CFD simulation
r=10; % rate at which the data is extracted
dt_r=dt*r;
%Time=20000*dt:dt_r:50000*dt;% flow time T = 100.5 - 1000
Time=20010*dt:dt_r:100000*dt;% flow time T = 100.5 - 1000
K=length(Time);
K_training = ceil(K/2); % Training dataset length
K_validation = K - K_training; % Validation dataset length
type=1;
%% WHICH DATASET TO LOAD - pressure or shear stress
  
%     a=a ;          % Temporal coefficients
%     phi=phi ;      % Eigenvectors (only ROM modes)
%     Y=Y ;      % Snapshot matrix of fluctuations
%     Y_mean=Y_mean; % mean velocity

    
    
%% CONTROL STUFF AND SIGNAL

fny=1/(2*dt_r); % juans was 5              %Frequence de Nyquist 1/(2*dt)
fphys= 0.15;   % found from FFT            %Frequence du ph�nomene physique
finput=fphys/fny;                   
u=idinput([K 4],'rbs', [ 0 finput]);       % signal to be used for excitation
%u=u*0.1;
%% H matrix

% m_H=1:50;
% H=a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)');
% Ident_check=round(a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)'),5);
% if Ident_check ~= eye(size(H))
%     msg = 'Identy matrix not obtained when A*inv(A) is calculated';
%     error(msg)  
% end 

%% data preparation
m=1:50;
%m=1:25; % number of POD modes to estimate % number of columns for transformation matirx
%H=H(m,:);
a=a(m,:); %Temporal coefficients  // there are 1800 coefficients in total
a=a(m,:);

% a is the same as a


%% System Identification w/ excitation

data = iddata(a(m,:)',u,dt_r,'Tstart',Time(1)); % create iddata Object for N4SID  // so here [] represents u (our input)
%n4sid

sysID_operation = input('Have you done sysID already? (yes=1, no=0): ');

if sysID_operation == 0 % Not computed

        
        Nx = input('Model order: '); % model order
        opt = n4sidOptions('InitialState','estimate','N4Weight','auto','Focus','Simulation','EstimateCovariance',true);
        [sysID,x0] = n4sid(data(1:K_training),Nx,'DisturbanceModel','estimate','CovarianceMatrix','estimate',opt); % stochastic identification

        if input('save system identification?')==1
            save([pathSID '9000_v3SYSID_Nx' num2str(Nx) '_M_' num2str( m(end)) '.mat' ],'sysID','x0') %saves system identification into control folder
        end 
        
elseif sysID_operation == 1 % Computed    
    
        model_order = input('What system order you want to import? ');
        %modes_input = input('what m to import?');
    
        load([pathSID '9000_v3SYSID_Nx' num2str(model_order) '_M_' num2str(m(end)) '.mat'] )
        disp(['SYSID_Nx' num2str(model_order) '_M_' num2str(m(end)) '.mat was loaded'])
else
    error('You did not input 1 or 0');
    return
    
    
end

A = sysID.A; % dynamics matrix
B = sysID.B; % input matrix
C = sysID.C; %output matrix

Obsv = obsv(A,C);
rank_Obsv = rank(Obsv);
if rank_Obsv == size(A,1)
    disp('The system is observable.');
else
    disp('The system is not observable.');
end

%%
% Plant matrices 
A = sysID.A; % dynamics matrix
B = sysID.B; % input matrix
C = sysID.C; % output matrix  
[Resid, AuCor]=resid(data(1:K_training),sysID);
e=Resid.y'; 
w=sysID.K*e;
%Q=cov(w',1); 
Q=sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
%%

% Qflavio = sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
% Q=Qflavio;
Nx = length(x0); % model order
disp(['max Q: ' num2str(max(max(Q)))]);

%% sensor_point are the indixes of the sensors found from QR pivoting
sensor_point = [15442 15420 16171 16230 4327 6449 20620 3674 22157 18987 7086 17617 3778 10490 4334 7803 21967 18203 3300 5567];
x_plotting2 = [x_plotting;x_plotting];
y_plotting2 = [y_plotting;y_plotting];
locations = [x_plotting2(sensor_point) y_plotting2(sensor_point)];
%sensor_point = [15420,15362,16230,16091,17294,18585,4093,19955,21032,19821];
p = length(sensor_point);
s = zeros(p,K);
s = Y(sensor_point,:); 

%% kalman filter
% S (must use the fish topology)
S = zeros(p,length(m));
S = phi(sensor_point,m); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
G = S*C;
% Measurements noise (Training dataset only)

%% noise measurement

measurement_noise = s(:,1:K_training) - S*a(m,1:K_training);%<<<<<<<<<<<<<<<<
R = cov(measurement_noise'); % covariance matrix
disp(['max R: ' num2str(max(max(R)))])


%plant = ss(A,[B eye(Nx)],S*C,[],dt_r);
plant = ss(A,[B eye(Nx)],S*C,[],dt_r);
N = zeros(Nx,p); % cross-covariance matrix
[kalmf,L,P] = kalman(plant,Q,R,N);

plant2 = ss(A,B,C,[],dt_r);

%% Nx = 100;
%% opt = n4sidOptions('InitialState','estimate','N4Weight','auto','Focus','Simulation','EstimateCovariance',true);
%% [sysIDv2,x0v2] = n4sid(datav2(1:900),Nx,'DisturbanceModel','estimate','CovarianceMatrix','estimate',opt); % stochastic identification

%% Control Design
% Q and R used for LQR is different than from Kalman
% matrix Q and R are positive definite and positive semi definite matrix, respectively
% Q_lqr = eye(Nx);
% want to minimise the energy of POD modes 
Q_lqr = transpose(C)*C;
R_lqr = 30*eye(4);
% R_lqr2(1,1) = 85;
% R_lqr2(2,2) = 90;
% R_lqr2(3,3) = 75;
% R_lqr2(4,4) = 65;

[K_lqr, S_lqr, e_lqr] = lqr(plant2,Q_lqr,R_lqr);

%% WHAT TO TAKE INTO STAR CCM

% A, B, C matrices from n4sid

% L matrix from Kalman filter

% K matrix from LQR

% S matrix 

% probably take in x_0 and a_e where these are the intialisation

save_matrix = input('Save matrices into csv files? (yes=1, no=0): ');

if save_matrix == 1 % Not computed
    writematrix(A, 'A2_9000_2.csv');
    writematrix(B, 'B2_9000_2.csv');
    writematrix(C, 'C2_9000_2.csv');
    writematrix(G, 'G2_9000_2.csv');
    writematrix(L, 'L2_9000_2.csv');
    writematrix(K_lqr, 'K2_9000_2.csv');
    writematrix(x0, 'x02_9000_2.csv');
elseif save_matrix == 0
    
else
end

% K_lqr3 = readtable('K2_9000_v3.csv');
% K_lqr3 = table2array(K_lqr3);
% 
% x_e = x0; % state initialisation
% a_e = C*x_e; % pre-allocation
% for k = 1:K
%     x_e = A*x_e + L*(s(:,k) - S*C*x_e);% - B*K_lqr3*x_e; % Remember this S includes phi*H
%     a_e(:,k) = C*x_e;
% end
% pathPOD='/home/wh219/FYP/Controller/Controller Save/';
% load([pathPOD 'POD_Control_v3_9000']);
% 
% n =11259;
% m_interest = 1:50;
% 
% x_length=[-3 11];
% y_length=[-4 4]; % Domain coordinates for sensor grid in star-ccm
%             
% for j = 1:16
%             snapshot=7660 + (j*10);
%             Magnitude_DNS = sqrt((Y(1:n,snapshot) + Y_mean(1:n)).^2 + (Y(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % DNS
% % 
% %             %SNAPSHOT TO COMPARE
% % 
% % 
% %             % plot
%                 Y_RECONSTRUCTED = zeros(size(Y_mean)); % pre-allocation
%                 for i = 1:length(m_interest)
%                     Y_RECONSTRUCTED = Y_RECONSTRUCTED + phi(:,m_interest(i))*a_e(m_interest(i),:); 
%                 end
%                 Magnitude_ROM = sqrt((Y_RECONSTRUCTED(1:n,snapshot) + Y_mean(1:n)).^2 + (Y_RECONSTRUCTED(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % ROM
% % 
% %                 DNS_contours = subplot(2,1,1);
% %                 zi = griddata(x_plotting,y_plotting,Magnitude_DNS,xq,yq);
% %                 contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
% %                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
% %                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
% %                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
% %                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
% %                 shading interp; axis equal; view(2); 
% %                 colormap(DNS_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
% %                 xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
% %                 xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
% %                 title(['DNS (timestep=' num2str(snapshot) ')'],'fontsize',14,'interpreter','latex');
% %                 set(gca,'linewidth',1);
% %                 set(gca,'TickLength',[0 0]);
% %                 set(gca,'Fontsize',12); 
% %                 %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
% % 
% %                 ESTIMATED_contours = subplot(2,1,2); 
% %                 zi = griddata(x_plotting,y_plotting,Magnitude_ROM,xq,yq);
% %                 contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
% %                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
% %                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
% %                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
% %                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
% %                 shading interp; axis equal; view(2); 
% %                 colormap(ESTIMATED_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
% %                 s=scatter(x_plotting(index_u2),y_plotting(index_u2),'filled','rs');  %% u velocity index being plotted for sensors
% %                 hold on;
% %                 s=scatter(x_plotting(index_v2),y_plotting(index_v2),'filled','bo');    %% v velocity index being plotted for sensors
% %                 hold off
% %                 xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
% %                 xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
% %                 title(['Reconstruction'],'fontsize',14,'interpreter','latex');
% %                 set(gca,'linewidth',1);
% %                 set(gca,'TickLength',[0 0]);
% %                 set(gca,'Fontsize',12); 
% %                 %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
% % 
%                 FITTING = subplot(4,4,j);
%                 FIT = 100*(1 - abs(Magnitude_ROM - Magnitude_DNS)./abs(Magnitude_DNS));
%                 zi = griddata(x_plotting,y_plotting,FIT,xq,yq);
%                 zi = max(zi,0); % remove all negative fit values
%                 contourf(xq,yq,zi,80,'Edgecolor','flat'); hold on;
%                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
%                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
%                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
%                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
%                 shading interp; axis equal; view(2); colormap(FITTING,jet); c = colorbar('eastoutside'); 
%                 caxis([0 100]); colorbar('Ticks',[0:25:100]);
%                 xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
%                 xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
%                 title(['Time = ',num2str(snapshot*0.1)],'fontsize',14,'interpreter','latex');
%                 set(gca,'linewidth',1);
%                 set(gca,'TickLength',[0 0]);
%                 set(gca,'Fontsize',12); 
%                 %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
%                 
% end               