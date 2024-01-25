%% Script for System Identification and Kalman v

clc; close all;
if input('Clear all variables?')==1
    clear; 
    path='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\';
    path2='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\Save Test\';
    pathSID='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\SID\';
    pathKalman='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\Kalman\';
    load([path 'POD_100000.mat'],'a','phi','Y_mean','Y')
    load([path2 'Griddata_100000.mat'],'xq','yq','x_plotting','y_plotting')
    %load([path 'fishsvd_200500_002TR40.mat']);

else
    clearvars -except X a phi Y_mean Y xq yq x_plotting y_plotting a_WS a_P Phi_P Phi_WS Y_P_mean Y_WS_mean Y_P Y_WS path pathSID pathKalman
 end


%% TEMPORAL DEFINITIONS
dt=0.01; %temporal resolution of CFD simulation
r=50; % rate at which the data is extracted
dt_r=dt*r;
%Time=20000*dt:dt_r:50000*dt;% flow time T = 100.5 - 1000
Time=10050*dt:dt_r:100000*dt;% flow time T = 100.5 - 1000
K=length(Time);
K_training = ceil(K/2); % Training dataset length
K_validation = K - K_training; % Validation dataset length
type=1;
%% WHICH DATASET TO LOAD - pressure or shear stress

workspace= input('Which dataset do you want to work with? (1=Shear Stress or 2=Pressure)');

if workspace == 1
%     a_fish=a_WS ;
%     Phi_fish=Phi_WS ;
%     Y_mat_fish=Y_WS ;
%     Y_mean_fish=Y_WS_mean;
%     plotname='Shear Stress';
    
    a_fish=a ;          % Temporal coefficients
    Phi_fish=phi ;      % Eigenvectors (only ROM modes)
    Y_mat_fish=Y ;      % Snapshot matrix of fluctuations
    Y_mean_fish=Y_mean; % mean velocity
    
%     
% elseif workspace == 2
%     a_fish=a_P ;
%     Phi_fish=Phi_P;
%     Y_mat_fish=Y_P ;
%     Y_mean_fish=Y_P_mean;
%     plotname='Pressure';
%     
% else 
%     error('ERROR: wrong selection input');
end 

%% H matrix

m_H=1:input('how many modes to estimate/construct H? flavio used 50 modes'); %% this will determine the number of sensor points
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
a=a(m,:); %Temporal coefficients  // there are 1800 coefficients in total
a_fish=a_fish(m,:);

% a_fish is the same as a


%% System Identification

data = iddata(a(m,:)',[],dt_r,'Tstart',Time(1)); % create iddata Object for N4SID  // so here [] represents u (our input)
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
[Resid, AuCor]=resid(data(1:K_training),sysID);
e=Resid.y'; 
w=sysID.K*e;
Q=cov(w',1); 
%Q=sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
%% Q matrix 
Q_=[];
for i=1:length(e)
    
    Q_(:,:,i)=w(:,i)*w(:,i)';
    
end 
Q_=sum(Q_,3)/length(e);
%%

% Qflavio = sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
% Q=Qflavio;
Nx = length(x0); % model order
disp(['max Q: ' num2str(max(max(Q)))])
%%
% C=[C; zeros(m_H(end)-m(end),length(C))];
%% Sensor selection & measurements
% SELECT THE SENSOR(S) TO INPUT INTO THE KALMAN PLANT
%sensor_point=100; % SELECT A NUMER, OR SET OF ROWS (SENSORS TO BE SELECTED)
sensor_point = [ceil(1:1:11259)]; % OR IF A SET OF SENSORS IS REQUIRED USE
sensor_point = [15420,15362,16230,16091,17294,18585,4093,19955,21032,19821];
load('QRsensor_points.mat' ,'QRsensors');
%sensor_point = QRsensors';
%VECTOR FORM []
% sensor points i think i need to change this to the index of the sensors
% found from the QR pivoiting algorithm
p = length(sensor_point);
s = zeros(p,K);
s = Y_mat_fish(sensor_point,:); 

%% kalman filter
% S (must use the fish topology)
S = zeros(p,length(m));
S = Phi_fish(sensor_point,m)*H(m,:); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Measurements noise (Training dataset only)

%% noise measurement

measurement_noise = s(:,1:K_training) - S*a(m,1:K_training);% g[k] in flavio's paper s[k] = S*a[k] + g[k]
R = cov(measurement_noise'); % covariance matrix
disp(['max R: ' num2str(max(max(R)))])


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
plant = ss(A,[B eye(Nx)],S*C,[],dt_r);
N = zeros(Nx,p); % cross-covariance matrix
[kalmf,L,P] = kalman(plant,Q,R,N);


%% Reconstruction of POD coefficients
% Estimation
x_e = x0; % state initialisation
a_e = C*x_e; % pre-allocation
for k = 1:K
    x_e = A*x_e + L*(s(:,k) - S*C*x_e); % Remember this S includes Phi_fish*H
    a_e(:,k) = C*x_e;
end

%only SID
x_e = x0; % state initialisation
a_SID = C*x_e; % pre-allocation

for k = 1:K
    x_e = A*x_e ; 
    a_SID(:,k) = C*x_e;
end


% FIT [%] modes
FIT_a_training = zeros(length(m),1); FIT_a_validation = FIT_a_training; % pre-allocation
for i = 1:length(m)
    % Training dataset
    FIT_a_training(i) = 100*(1-((norm(a(i,1:K_training) - a_e(i,1:K_training)))/(norm(a(i,1:K_training) - mean(a(i,1:K_training))))));
    
    % Validation dataset
    FIT_a_validation(i) = 100*(1-((norm(a(i,K_training+1:end) - a_e(i,K_training+1:end)))/(norm(a(i,K_training+1:end) - mean(a(i,K_training+1:end))))));
end
FIT_a = [FIT_a_training, FIT_a_validation];  

%% Reconstruction of u' and v' fields

% run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\FIT_U_V.m')
% % 
% % %% Reconstruction video
% run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\Reconstruction.m')
% %  %% POD COEFFICIENTS PLOT1
% run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\POD_coeff_plots.m')

% %% Save results
% if input('save result')==1
%     name_save = input('Save reconstruction? (if so, type name): ','s');
%     save([pathSID name_save],'Nx','p','sensor_point','a_e','FIT_a','FIT_u_training','FIT_u_validation','FIT_v_training','FIT_v_validation','Y_ROM');
% %     save([pathSID name_save],'Nx','p','sensor_point','a_e','FIT_a');
% end 

type=2;
