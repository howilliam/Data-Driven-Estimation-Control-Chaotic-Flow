%% Iterative search code
clc; close all;
if input('Clear all variables?')==1
    clear; 
    path='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\';
    pathSID='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\SID\';
    pathKalman='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\Kalman\';
    load([path 'POD_HPC_200500_002TR20.mat'],'a','phi','Y_mean','Y');
    load([path 'Griddata_HPC_200500_002TR20.mat'],'xq','yq','x_plotting','y_plotting')
    load([path 'fishsvd_200500_002TR40.mat']);

else
    clearvars -except X a phi Y_mean Y xq yq x_plotting y_plotting a_WS a_P Phi_P Phi_WS Y_P_mean Y_WS_mean Y_P Y_WS path pathSID pathKalman
 end


%% TEMPORAL DEFINITIONS
dt=0.002; %temporal resolution of CFD simulation
r=40; % rate at which the data is extracted
dt_r=dt*r;
Time=100020*dt:dt_r:250220*dt;% flow time

K=length(Time);
K_training = ceil(K*0.7); % Training dataset length
K_validation = K - K_training; % Validation dataset length
type=1;
%% WHICH DATASET TO LOAD - pressure or shear stress

workspace= input('Which dataset do you want to work with? (1=Shear Stress or 2=Pressure)');

if workspace == 1
    a_fish=a_WS ;
    Phi_fish=Phi_WS ;
    Y_mat_fish=Y_WS ;
    Y_mean_fish=Y_WS_mean;
    plotname='Shear Stress';
    
elseif workspace == 2
    a_fish=a_P ;
    Phi_fish=Phi_P;
    Y_mat_fish=Y_P ;
    Y_mean_fish=Y_P_mean;
    plotname='Pressure';
    
else 
    error('ERROR: wrong selection input');
end 

%% H matrix

m_H=1:input('how many modes to estimate/construct H?');
H=a_fish(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)');
Ident_check=round(a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)'),5);
if Ident_check ~= eye(size(H))
    msg = 'Identy matrix not obtained when A*inv(A) is calculated';
    error(msg)  
end 

%% data preparation
m=m_H;
%m=1:25; % number of POD modes to estimate % number of columns for transformation matirx
H=H(m,:);
a=a(m,:);
a_fish=a_fish(m,:);



%% System Identification

data = iddata(a(m,:)',[],dt_r,'Tstart',Time(1)); % create iddata Object for N4SID
%n4sid

sysID_operation = input('Have you done sysID already? (yes=1, no=0): ');

if sysID_operation == 0 % Not computed

        
        Nx = input('Model order: '); % model order
        opt = n4sidOptions('InitialState','estimate','N4Weight','auto','Focus','Simulation','EstimateCovariance',true);
        [sysID,x0] = n4sid(data(1:K_training),Nx,'DisturbanceModel','estimate','CovarianceMatrix','estimate',opt); % stochastic identification

        if input('save system identification?')==1
            save([pathSID 'SYSID_Nx' num2str(Nx) '_M_' num2str( m(end)) '.mat' ],'sysID','x0')
        end 
        
elseif sysID_operation == 1 % Computed    
    
        model_order = input('What system order you want to import? ');
        %modes_input = input('what m to import?');
    
        load([pathSID 'SYSID_Nx' num2str(model_order) '_M_' num2str(m(end)) '.mat'] )
        disp(['SYSID_Nx' num2str(model_order) '_M_' num2str(m(end)) '.mat was loaded'])
else
    error('You did not input 1 or 0');
    return
    
    
end

%%
% Plant matrices 
A = sysID.A; % dynamics matrix
B = sysID.B; % input matrix
C = sysID.C; % output matrix
Nx = length(x0); % model order


% SOLUTION 1
[Resid, AuCor]=resid(data(1:K_training),sysID);
e=Resid.y'; 
w=sysID.K*e;
Q=cov(w',1); 

% SOLUTION 2
Q=sysID.K*sysID.NoiseVariance*sysID.K';
%% Sensor selection & measurements
%sensor_point = ceil((500-1).*rand([50,1]) + 1);% indecies of sensors selected from snapshot matrix 
avgs=[];
goodsens=[];
for J=1:500
m=m_H;
sensor_point = [J];


p = length(sensor_point);
s = zeros(p,K);
s = Y_mat_fish(sensor_point,:); 

% kalman filter
% S (must use the fish topology)
S = zeros(p,length(m));
S = Phi_fish(sensor_point,m)*H; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Measurements noise (Training dataset only)

% noise measurement

measurement_noise = s(:,1:K_training) - S*a(m,1:K_training);%<<<<<<<<<<<<<<<<
R = cov(measurement_noise'); % covariance matrix
%Q=0.00001*Q;
% % R=0.00000001*R;
% Cross-covariance (Training dataset only)
% random process noise vector 
% K_LD= length(A);
% w = randn(Nx,K_LD); 
% w = chol(Q,'lower')*w; 
% 
% for i = 1:Nx
%     for j = 1:length(sensor_point) 
%         N(i,j) = (1/K_LD)*w(i,:)*measurement_noise(j,:)';
%     end
% end

% Kalman filter gain L
plant = ss(A,[B eye(Nx)],S*C,[],dt);
N = zeros(Nx,p); % cross-covariance matrix
[kalmf,L,P] = kalman(plant,Q,R,N);

%     A = kalmf.A; % dynamics matrix
%     B = kalmf.B; % input matrix
%     C = kalmf.C; % output matrix

% Reconstruction of POD coefficients
% Estimation
x_e = x0; % state initialisation
a_e = C*x_e; % pre-allocation
for k = 1:K
    x_e = A*x_e + L*(s(:,k) - S*C*x_e);
    a_e(:,k) = C*x_e;
end

% FIT [%] modes
m_int=1:5;
FIT_a_training = zeros(length(m_int),1); FIT_a_validation = FIT_a_training; % pre-allocation

for i = m_int
    % Training dataset
    FIT_a_training(i,:) = 100*(1-((norm(a(i,1:K_training) - a_e(i,1:K_training)))/(norm(a(i,1:K_training) - mean(a(i,1:K_training))))));
    
    % Validation dataset
    FIT_a_validation(i,:) = 100*(1-((norm(a(i,K_training+1:end) - a_e(i,K_training+1:end)))/(norm(a(i,K_training+1:end) - mean(a(i,K_training+1:end))))));
end
FIT_a = [FIT_a_training, FIT_a_validation];  
FIT_a_all(1,:,J) = FIT_a_training;FIT_a_all(2,:,J) = FIT_a_validation;
%average search
avgs(J,:) = [J mean(FIT_a(m_int,1)) mean(FIT_a(m_int,2))]; 

%threshold search
THRESHOLD=87;
if sum(FIT_a(1:m_int(end),1)>THRESHOLD)== length(1:5)
    goodsens=[goodsens; J ]
end

end %%end for loop J , sensor search

%threshold sensors
goodsens
sensorplot(goodsens,['Nx 20, ' num2str(THRESHOLD) '\% FIT threshhold,' plotname ' sensors'])
% best sensors in VALIDATION DATASET
sorted_avgs=sortrows(avgs,3,'descend');
TOP_20=sorted_avgs(1:20,:)

%%
% if input('Save results?')==1
% save([pathSID 'Iterative_Nx' num2str(Nx) '_M_' num2str(m(end)) '.mat'], 'avgs', 'TOP_20', 'Nx' )
% end

% if input('plot sensor location?')==1
%     run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\supporting\Plotsens.m')
% end
%%
figure
tiledlayout(2,1)
nexttile
%subplot(2,1,1)%SENSOR LOCATION
    load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\fish_sorted_HPC_200500_002TR20.mat','Y_fish','X_fish');% contains data order as used in codes
    plot(X_fish([1:2:500 500]),Y_fish([1:2:500 500]),'k', 'linewidth', 1.3); hold on;
    plot(X_fish([1 2:2:500]),Y_fish([1 2:2:500]),'k', 'linewidth', 1.3); axis equal
    ylim([9.5 10.5]); 
    xlim([8 12])
    set(gca,'XTick',[], 'YTick', [])
    title(['Fish, Top 20 ' plotname ' sensors (by average)'],'interpreter','latex','fontsize',15)
    s=scatter(X_fish(TOP_20(:,1)),Y_fish(TOP_20(:,1)),'ro')

%subplot(2,1,2) %BAR GRAPH
nexttile
    
    bar((X_fish(1:1:end)-X_fish(1))/(X_fish(end)-X_fish(1)),avgs(:,3),3000); hold on;
    % plot(X_fish(2:2:end)/X_fish(end),avgs(2:2:end,3),'k','linewidth',1.3)
    ylim([0 100])
    xlabel('x/L','interpreter','latex')
    ylabel('Average FIT [\%]','interpreter','latex');
    title(['Average FIT of modes:  ' plotname ' sensors'],'interpreter','latex','fontsize',15)

%save plot
if input('Save plot?')==1
    f=gcf;
    exportgraphics(f,[path 'images\iterative\ORDER ' num2str(Nx) ' M ' num2str(m(end)) ' ITERATIVE - ' plotname '.jpg'],'Resolution',300)
    exportgraphics(f,[path 'images\iterative\ORDER' num2str(Nx) '_M' num2str(m(end)) '_ITERATIVE_' plotname '.eps'],'ContentType','vector')
end