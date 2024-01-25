%% Script for System Identification and Kalman v

clc; close all;
if input('Clear all variables?')==1
    clear; 
    path='/home/wh219/FYP/Open Loop/'; %path to POD 
    path2='/home/wh219/FYP/Open Loop/';          % path to grid data                open loop scenario
    pathSID='/home/wh219/FYP/Open Loop/'; %path to save system identifification
    pathKalman='/home/wh219/FYP/Open Loop/'; %path to save Kalman
    
%     path='/home/wh219/FYP/Controller/Controller Save/'; %path to POD 
%     path2='/home/wh219/FYP/Controller/Controller Save/';          % path to grid data
%     pathSID='/home/wh219/FYP/Controller/Controller Save/'; %path to save system identifification
%     pathKalman='/home/wh219/FYP/Controller/Controller Save/'; %path to save Kalman
    
    load([path 'POD_Open_9000'],'a','phi','Y_mean','Y')
    load([path2 'Griddata_Open_9000.mat'],'xq','yq','x_plotting','y_plotting')
    %load([path 'fishsvd_200500_002TR40.mat']);

else
    clearvars -except X a phi Y_mean Y xq yq x_plotting y_plotting a_WS a_P Phi_P Phi_WS Y_P_mean Y_WS_mean Y_P Y_WS path pathSID pathKalman
 end


%% TEMPORAL DEFINITIONS
dt=0.01; %temporal resolution of CFD simulation
r=10; % rate at which the data is extracted
dt_r=dt*r;
%Time=20000*dt:dt_r:50000*dt;% flow time T = 100.5 - 1000
Time=10010*dt:dt_r:100000*dt;% flow time T = 100.5 - 1000
K=length(Time);
K_training = ceil(K/2); % Training dataset length
K_validation = K - K_training; % Validation dataset length
n = length(x_plotting);
type=1;
%% WHICH DATASET TO LOAD - pressure or shear stress

workspace= input('Which dataset do you want to work with? (1=Shear Stress or 2=Pressure)');

if workspace == 1
%     a=a_WS ;
%     phi=Phi_WS ;
%     Y_mat_fish=Y_WS ;
%     Y_mean=Y_WS_mean;
%     plotname='Shear Stress';
    
%     a=a ;          % Temporal coefficients
%     phi=phi ;      % Eigenvectors (only ROM modes)
%     Y=Y ;      % Snapshot matrix of fluctuations
%     Y_mean=Y_mean; % mean velocity
    
%     
% elseif workspace == 2
%     a=a_P ;
%     phi=Phi_P;
%     Y_mat_fish=Y_P ;
%     Y_mean=Y_P_mean;
%     plotname='Pressure';
%     
% else 
%     error('ERROR: wrong selection input');
end 

%% H matrix

% m_H=1:input('how many modes to estimate/construct H? flavio used 50 modes'); %% this will determine the number of sensor points
% H=a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)');
% Ident_check=round(a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)'),5);
% if Ident_check ~= eye(size(H))
%     msg = 'Identy matrix not obtained when A*inv(A) is calculated';
%     error(msg)  
% end 

%% data preparation
%m=m_H;
m=1:50; % number of POD modes to estimate % number of columns for transformation matirx
%H=H(m,:);
% a=a(m,:); %Temporal coefficients  // there are 1800 coefficients in total
% a=a(m,:);

% a is the same as a


%% System Identification

data = iddata(a(m,:)',[],dt_r,'Tstart',Time(1)); % create iddata Object for N4SID  // so here [] represents u (our input)
%n4sid

sysID_operation = input('Have you done sysID already? (yes=1, no=0): ');

if sysID_operation == 0 % Not computed

        
        Nx = input('Model order: '); % model order
        opt = n4sidOptions('InitialState','estimate','N4Weight','auto','Focus','Simulation','EstimateCovariance',true);
        [sysID,x0] = n4sid(data(1:K_training),Nx,'DisturbanceModel','estimate','CovarianceMatrix','estimate',opt); % stochastic identification

        if input('save system identification?')==1
            save([pathSID 'Excite_SYSID_Nx_' num2str(Nx) '_M_' num2str( m(end)) '.mat' ],'sysID','x0')
        end 
        
elseif sysID_operation == 1 % Computed    
    
        model_order = input('What system order you want to import? ');
        %modes_input = input('what m to import?');
    
        load([pathSID '9000_SYSID_Nx_' num2str(model_order) '_M_' num2str(m(end)) '.mat'] )
        disp(['SYSID_Nx' num2str(model_order) '_M_' num2str(m(end)) '.mat was loaded'])
else
    error('You did not input 1 or 0');
    return
    
    
end

%%
% Plant matrices 
A = sysID.A; % dynamics matrix
B = sysID.B; % input matrix
C = sysID.C; % output matrix  // this is my control law C, same that controls the plant
[Resid, AuCor]=resid(data(1:K_training),sysID);
e=Resid.y'; 
w=sysID.K*e;
Q=cov(w',1); 
%Q=sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
%%
% Qflavio = sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
% Q=Qflavio;
Nx = length(x0); % model order
disp(['max Q: ' num2str(max(max(Q)))])
%% this section previously used sensor selction that was random, go straight to QR to pick sensor location
% %%
% % C=[C; zeros(m_H(end)-m(end),length(C))];
% %% Sensor selection & measurements
% % SELECT THE SENSOR(S) TO INPUT INTO THE KALMAN PLANT
% %sensor_point=100; % SELECT A NUMER, OR SET OF ROWS (SENSORS TO BE SELECTED)
% sensor_point = [ceil(1:1:100)]; % OR IF A SET OF SENSORS IS REQUIRED USE
% sensor_point = [15420,15362,16230,16091,17294,18585,4093,19955,21032,19821];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%load('QRsensor_points.mat' ,'QRsensors');
% %sensor_point = QRsensors';
% %VECTOR FORM []
% % sensor points i think i need to change this to the index of the sensors
% % found from the QR pivoiting algorithm
% p = length(sensor_point);
% s = zeros(p,K);
% s = Y_mat_fish(sensor_point,:); 
% 
% %% kalman filter
% % S (must use the fish topology)
% S = zeros(p,length(m));
% S = phi(sensor_point,m); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% % Measurements noise (Training dataset only)
% 
% %% noise measurement
% 
% measurement_noise = s(:,1:K_training) - S*a(m,1:K_training);% g[k] in flavio's paper s[k] = S*a[k] + g[k]
% R = cov(measurement_noise'); % covariance matrix
% disp(['max R: ' num2str(max(max(R)))])


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
% plant = ss(A,[B eye(Nx)],S*C,[],dt_r);
% N = zeros(Nx,p); % cross-covariance matrix
% [kalmf,L,P] = kalman(plant,Q,R,N);


%% Build the Kalman filter by selecting sensors from QR pivoting 
type=2;

%% QR pivoting 

%Plant is loaded before

%number of POD MODES
p_modes=length(m);

%% QR pivoting

    % DT Lyapunov equation for the controllability Gramian 
    Wc = dlyap(sysID.A,sysID.K*(sysID.K)');
    
    % Output controllability Gramian
    Woc = sysID.C*Wc*(sysID.C)';
    %Woc = sysID.C*Wc*(sysID.C)';
    
        % Check positive-definite property
    if norm(Woc - Woc') < 1e-9 && min(real(eig(Woc))) > -1e-9
        fprintf(1,'The output controllability Gramian is positive-definite \n');
    end
%%    
    % Cholesky decomposition: Woc = L*L'
    [Lc,flag] = chol(Woc,'lower'); 
    if flag == 0 && norm(Lc*Lc' - Woc) < 1e-9
        fprintf(1,"The lower triangular factor satisfies L*L' - Woc = 0, within roundoff error %e \n",norm(Lc*Lc' - Woc));
    else 
        fprintf(1,'Error in Cholesky decomposition \n');
    end 
    
    % QR pivoting algorithm
    [~,R_qr,pivot] = qr((phi(:,1:p_modes)*Lc)','vector'); 
    
    % Save optimally ranked sensors
     %QRsensors = pivot(1:p_modes)'; %will have same number of sensors as modes
    QRsensors = pivot(1:20)'; %will have same number of sensors as modes
    
    %some index exceeds 11259 becuase it looks at u and v in the combined
    %matrix
    
    storage = QRsensors > length(x_plotting);
    
    QRsensors_v = [];
    QRsensors_u = [];
    
    for i = 1:length(QRsensors)
        if storage(i) == 1
            QRsensors_v(end+1) = QRsensors(i);
        else
            QRsensors_u(end+1) = QRsensors(i);
        end
    end 
            
    save(['QRsensors_',num2str(Nx)],'QRsensors');
    %sensorplot2(QRsensors,['Nx=' num2str(Nx) ', QR best ranked ' plotname ' sensors'],p_modes);
    sensorplot2(QRsensors_u,QRsensors_v,['Nx=' num2str(Nx) ', QR best ranked sensors'],p_modes);
    %%
    %Oversensed system
    %[~,~,pivot2] = qr((Phi_fish(:,1:p_modes)*Lc)*(Phi_fish(:,1:p_modes)*Lc)','vector');
    
    
    %sensorplot(pivot2(1:20),['Nx=' num2str(Nx) ', oversensed-QR best ranked ' plotname ' sensors'],20);
    
    %S_optimal = S(sort(QRsensors,'ascend'),:);
    S_optimal = phi(sort(QRsensors,'ascend'),m);
    % Measurement noise (Training dataset only)
    % s_optimal = s(sort(QRsensors,'ascend'),:); % optimal sensor inputs 
    s_optimal = Y(sort(QRsensors,'ascend'),:); % optimal sensor inputs 
    measurement_noise = s_optimal(:,1:K_training) - S_optimal*a(m,1:K_training);
    R_optimal = cov(measurement_noise'); % covariance matrix
    
    
    
    % Kalman filter gain L
plant_optimal = ss(A,[B eye(Nx)],S_optimal*C,[],dt_r);
N_optimal = zeros(Nx,length(QRsensors)); % noise cross covariance
[kalmf_optimal,L_optimal,P_optimal] = kalman(plant_optimal,Q,R_optimal,N_optimal);

%% Estimated POD coefficients from optimal sensors 
x_e = x0; % state initialisation
a_e = C*x_e; % pre-allocation
for k = 1:K
    x_e = A*x_e + L_optimal*(s_optimal(:,k) - S_optimal*C*x_e);
    a_e(:,k) = C*x_e;
end



% FIT [%]
FIT_a_training = zeros(length(m),1); FIT_a_validation = FIT_a_training;
for i = 1:length(m)
    % Training dataset
    FIT_a_training(i) = 100*(1-((norm(a(i,1:K_training) - a_e(i,1:K_training)))/(norm(a(i,1:K_training) - mean(a(i,1:K_training))))));
    
    % Validation dataset
    FIT_a_validation(i) = 100*(1-((norm(a(i,K_training+1:end) - a_e(i,K_training+1:end)))/(norm(a(i,K_training+1:end) - mean(a(i,K_training+1:end))))));
end
FIT_a_optimal = [FIT_a_training, FIT_a_validation];

index_u2 = QRsensors_u;    
index_v2 = QRsensors_v - n;   
% Plot POD coefficients reconstruction
run('/home/wh219/FYP/Codes/MSCproj/200500 data/supporting/FIT_U_V.m')
run('/home/wh219/FYP/Codes/MSCproj/200500 data/supporting/Reconstructionv3.m')
run('/home/wh219/FYP/Codes/MSCproj/200500 data/supporting/POD_coeff_plots.m')
