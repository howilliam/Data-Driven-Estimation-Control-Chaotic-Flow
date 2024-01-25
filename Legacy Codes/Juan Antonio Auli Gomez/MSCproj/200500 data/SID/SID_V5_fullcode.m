%% kalmandrafts
clear; 
clc; close all;

path='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\';
path2='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\Save Test\';
pathSID='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\SID\';
pathKalman='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\Kalman\';
load([path 'POD_100000.mat'],'a','phi','Y_mean','Y')
load([path2 'Griddata_100000.mat'],'xq','yq','x_plotting','y_plotting')

% path='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\';
% pathSID='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\SID\';
% pathKalman='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\Kalman\';
% load([path 'POD_HPC_400500_002TR20.mat'],'a','phi','Y_mean','Y');
% load([path 'Griddata_HPC_400500_002TR20.mat'],'xq','yq','x_plotting','y_plotting')

% load([path 'fishsvd_400500_002TR40.mat']);

%% TEMPORAL DEFINITIONS
dt=0.01; %temporal resolution of CFD simulation
r=50; % rate at which the data is extracted
dt_r=dt*r;
Time=10050*dt:dt_r:100000*dt;% flow time
n = length(x_plotting);

K=length(Time);
K_training = ceil(K*0.7); % Training dataset length
K_validation = K - K_training; % Validation dataset length
type=1;
%% WHICH DATASET TO LOAD - pressure or shear stress

% workspace= input('Which dataset do you want to work with? (1=Shear Stress or 2=Pressure)');
% 
% if workspace == 1
%     a=a_WS ;
%     phi=Phi_WS ;
%     Y=Y_WS ;
%     Y_mean=Y_WS_mean;
%     plotname='Shear Stress';
    
% elseif workspace == 2
%     a=a_P ;
%     phi=Phi_P;
%     Y=Y_P ;
%     Y_mean=Y_P_mean;
%     plotname='Pressure';
%     
% else 
%     error('ERROR: wrong selection input');
% end 

%% H matrix

% m_H=1:input('how many modes to estimate/construct H?');
% H=a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)');
% Ident_check=round(a(m_H,1:K_validation)*a(m_H,1:K_validation)'*inv(a(m_H,1:K_validation)*a(m_H,1:K_validation)'),5);
% if Ident_check ~= eye(size(H))
%     msg = 'Identy matrix not obtained when A*inv(A) is calculated';
%     error(msg)  
% end 
% 
% %% data preparation
% m=m_H;
% %m=1:25; % number of POD modes to estimate % number of columns for transformation matirx
% H=H(m,:);
% a=a(m,:);
% a=a(m,:);



%% System Identification
m = 1:50; % number of POD modes to estimate

data = iddata(a(m,:)',[],dt_r,'Tstart',Time(1)); % create iddata Object for N4SID
%n4sid

sysID_operation = input('Have you done sysID already? (yes=1, no=0): ');

if sysID_operation == 0 % Not computed

        data = iddata(a(m,:)',[],dt_r,'Tstart',Time(1)); % create iddata Object for N4SID
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
% [Resid, AuCor]=resid(data(K_training+1:end),sysID);
% e=Resid.y'; 
% w=sysID.K*e;
% Q=cov(w',1); 
Q=sysID.K*sysID.NoiseVariance*sysID.K';
%% Q matrix 
% Q_=[];
% for i=1:length(e)
%     
%     Q_(:,:,i)=w(:,i)*w(:,i)';
%     
% end 
% Q_=sum(Q_,3)/length(e);
%%

%Qflavio = sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance

Nx = length(x0); % model order
% disp(['max Q: ' num2str(max(max(Q)))])
%%
% C=[C; zeros(m_H(end)-m(end),length(C))];
%% Sensor selection & measurements
%sensor_point = ceil((500-1).*rand([50,1]) + 1)% indecies of sensors selected from snapshot matrix 
%sensor_point = [452;446;168;350;100;17;373;251;241;453;306;310;430;403;289;93;121;444;16;246;85;490;357;251;237;31;342;23;37;262;50;410;409;362;76;331;260;487;325;401;228;217;413;43;68;88;197;416;402;32]
%sensor_point=[289;342;274;214;323;325;340;319;473;106;355;119;61;305;226;230;332;386;176;332;209;422;417;129;308;292;271;436;134;160;61;470;324;241;321;273;325;273;361;262;497;111;54;56;33;203;225;184;382;315]
%sensor_point = [1,177,327,400,101,103,459,461,175,413,329,476,411,489,474,478,2,284,286,472];
sensor_point = [ceil(1:1:11259)];

p = length(sensor_point);
s = zeros(2*p,K);
s(1:2:end,:) = Y(sensor_point,:); % u'
s(2:2:end,:) = Y(sensor_point+n,:); % v'
% s = Y(sensor_point,:); 

%% kalman filter
% S (must use the fish topology)
S = zeros(2*p,length(m));
% S = phi(sensor_point,m)*H(m,:); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Measurements noise (Training dataset only)
S(1:2:end,:) = phi(sensor_point,m); % u'
S(2:2:end,:) = phi(sensor_point+n,m); % v'

%% noise measurement

measurement_noise = s(:,1:K_training) - S*a(m,1:K_training);%<<<<<<<<<<<<<<<<
R = cov(measurement_noise'); % covariance matrix

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
N = zeros(Nx,2*p); % cross-covariance matrix
[kalmf,L,P] = kalman(plant,Q,R,N);


%% Reconstruction of POD coefficients
% Estimation
% x_e = x0; % state initialisation
% a_e = C*x_e; % pre-allocation
% for k = 1:K
%     x_e = A*x_e + L*(s(:,k) - S*C*x_e); % Remember this S includes phi*H
%     a_e(:,k) = C*x_e;
% end
% 
% % FIT [%] modes
% FIT_a_training = zeros(length(m),1); FIT_a_validation = FIT_a_training; % pre-allocation
% for i = 1:length(m)
%     % Training dataset
%     FIT_a_training(i) = 100*(1-((norm(a(i,1:K_training) - a_e(i,1:K_training)))/(norm(a(i,1:K_training) - mean(a(i,1:K_training))))));
%     
%     % Validation dataset
%     FIT_a_validation(i) = 100*(1-((norm(a(i,K_training+1:end) - a_e(i,K_training+1:end)))/(norm(a(i,K_training+1:end) - mean(a(i,K_training+1:end))))));
% end
% FIT_a = [FIT_a_training, FIT_a_validation];  
% 
% %% Reconstruction of u' and v' fields
% 
% % run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\FIT_U_V.m')
% run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\FIT_U_V.m')
% %% Reconstruction video
% % run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\Reconstruction.m')
% run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\Reconstruction.m')
% %% POD COEFFICIENTS PLOT
% % run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\POD_coeff_plots.m')
% % run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\POD_coeff_plots.m')
% %% Save results
% if input('save result')==1
%     name_save = input('Save reconstruction? (if so, type name): ');
%     save([pathSID name_save],'Nx','p','sensor_point','a_e','FIT_a','FIT_u_training','FIT_u_validation','FIT_v_training','FIT_v_validation','Y_ROM'); 
% end 

type=2;
%% Sparse sensor placement 

%Balanced trucation 
% kalman filter Discr Time system
kalman_filter = ss(A-L*S*C,... %Kalman A matrix , remember that this S contains also thr H matrix
    L,...                      %Kalman B matrix     
    S*C,...                    %Kalman C matrix
    [],...                     %Kalman D matrix
    dt_r);                     %Temooral def

% Stability
figure('name','Stability of Kalman filter');
[eig_kalman_filter] = eig(kalman_filter.A);
theta = 0:0.01:2*pi; xc = cos(theta); yc = sin(theta);
plot(xc,yc,'k--','Linewidth',1); hold on; % unit circle
plot(real(eig_kalman_filter),imag(eig_kalman_filter),'r.','Markersize',20);
grid on; axis equal; xlabel('Re'); ylabel('Im'); 
xlim([-1.1 1.1]); ylim([-1.1 1.1]); xticks(-1:0.5:1); yticks(-1:0.5:1);
title('Eigenvalues of the Kalman filter');

% Unbalanced Gramians 
Wc = dlyap(kalman_filter.A,kalman_filter.B*kalman_filter.B');
Wo = dlyap(kalman_filter.A',kalman_filter.C'*kalman_filter.C);

[sysb,g,Ti,T] = balreal(kalman_filter); % balanced system

% Balanced Gramians (equal and diagonal)
Wc_b = dlyap(sysb.A,sysb.B*sysb.B');
Wo_b = dlyap(sysb.A',sysb.C'*sysb.C);

% Hankel singular values
figure('name','Balanced truncation'); 
subplot(1,2,1);
plot(g./sum(g),'r.-','Markersize',20,'Linewidth',2);
xlabel('State'); ylabel('Normalised HSVs'); grid on;
subplot(1,2,2);
plot(cumsum(g./sum(g)),'r.-','Markersize',20,'Linewidth',2);
xlabel('State'); ylabel('Cumulative Normalised HSVs'); grid on;

r = input('Sensor budget r: ');
PSI_r = Ti(:,1:r);
PHI_star_r = T(1:r,:);

% Balanced-truncated system
sys_bt = ss(PHI_star_r*kalman_filter.A*PSI_r,... % A matrix
    PHI_star_r*kalman_filter.B,... % B matrix
    kalman_filter.C*PSI_r,... % C matrix
    [],... % D matrix
    dt_r);

%% Optimal sensor placement
% QR pivoting
[Qqr,Rqr,Pqr] = qr(sys_bt.B);
   
% Selection matrix
S_B = Pqr(:,1:r); % optimal selection matrix
B_hat = kalman_filter.B*S_B; % new actuator matrix

% Find the indices of the optimal sensors
for i = 1:r
    for j = 1:p
        found(i,j) = norm(B_hat(:,i) - kalman_filter.B(:,j));
    end
    optimal_input_index(i) = find(found(i,:) == 0); % contains rows of optimal sensors
end

% figure('name','Optimal sensors');
% subplot(5,6,[1:5, 7:11]); % B matrix 
% h = imagesc(kalman_filter.B); colormap pink; hold on;
% set(gca,'ytick',[]); xticks(sort(optimal_input_index,'ascend')); xtickangle(90);
% ad = 0.2*ones(size(kalman_filter.B)); ad(:,optimal_input_index) = 1; set(h,'AlphaData',ad);
% xlabel('Number of $2q$ sensor inputs','fontsize',14); ylabel('Number of model states $N_{X}$','fontsize',14); 
% title('Optimal $\qhbf{B}$ matrix columns','fontsize',14);
% 
% subplot(5,6,[6:6:30]); % s input vector 
% h = imagesc(repmat(rand(length(sensor_point),1),[1 2])); axis equal; colormap pink; hold on;
% set(gca,'xtick',[]); yticks(sort(optimal_input_index,'ascend'));
% xlim([0.5 2.5]); ylim([0.5 length(sensor_point)+0.5]);
% ad = 0.1*ones(length(sensor_point),2); ad(optimal_input_index,:) = 1; set(h,'AlphaData',ad);
% ylabel('Number of $2q$ sensor inputs','fontsize',14); 
% title('Input $\mathbf{s}$','fontsize',14);



%% Best sensors
% Find optimal sensors where pressure or shear stress is extracted
optimal_sensor_location = ceil(sort(optimal_input_index,'ascend')/2);  %%%% not /2, fish sens
[~,w] = unique(optimal_sensor_location,'stable'); 
optimal_sensor_location_uv = optimal_sensor_location(setdiff(1:numel(optimal_sensor_location),w)); % u',v'
clear w;


index_sort = sort(optimal_input_index,'ascend'); 
for i = 1:length(index_sort)
    if mod(index_sort(i),2) == 0 && ismember(optimal_sensor_location(i), optimal_sensor_location_uv) == 0 % even index
        optimal_sensor_location_v(i) = optimal_sensor_location(i); % v' 
    elseif mod(index_sort(i),2) ~= 0 && ismember(optimal_sensor_location(i), optimal_sensor_location_uv) == 0 % odd index
        optimal_sensor_location_u(i) = optimal_sensor_location(i); % u'
    end
end
% storage = optimal_input_index > length(x_plotting);
%     
%     QRsensors_v = [];
%     QRsensors_u = [];
%     
%     for i = 1:length(optimal_input_index)
%         if storage(i) == 1
%             QRsensors_v(end+1) = optimal_input_index(i);
%         else
%             QRsensors_u(end+1) = optimal_input_index(i);
%         end
%     end 
   
%sensorplot(optimal_input_index,'optimal input index fish sensors')
%sensorplot2(,['Optimal location: ' num2str(r) ' sensors'],length(optimal_input_index))
%         saveas(gcf,[pathvid 'Optimal sensor' num2str(r) '_'  plotname],'epsc') %save images
%         saveas(gcf,[pathvid 'Optimal sensor' num2str(r) '_' plotname])


optimal_sensor_location_u = optimal_sensor_location_u(optimal_sensor_location_u~=0);
optimal_sensor_location_v = optimal_sensor_location_v(optimal_sensor_location_v~=0);
best_sensors_u = sensor_point(optimal_sensor_location_u); 
best_sensors_v = sensor_point(optimal_sensor_location_v); 
best_sensors_uv = sensor_point(optimal_sensor_location_uv); 
clear optimal_sensor_location optimal_sensor_location_u optimal_sensor_location_v optimal_sensor_location_uv index_sort; 

figure('name','Optimal sensor library');
p1 = plot(x_plotting,y_plotting','k.','Markersize',1); hold on; % mesh points
p2 = plot(x_plotting(sensor_point),y_plotting(sensor_point),'r.','Markersize',12); hold on; % original sensor library
p3 = plot(x_plotting(best_sensors_u),y_plotting(best_sensors_u),'rx','linewidth',3,'Markersize',10); hold on; % u' sensor
p4 = plot(x_plotting(best_sensors_v),y_plotting(best_sensors_v),'bx','linewidth',3,'Markersize',10); hold on; % v' sensor
p5 = plot(x_plotting(best_sensors_uv),y_plotting(best_sensors_uv),'kx','linewidth',3,'Markersize',10); % u',v' sensors
grid off; axis equal; view(2);
xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
legend([p1,p2,p3,p4,p5],'Mesh','Original sensor library','Optimal sensors $(u)$','Optimal sensors $(v)$','Optimal sensors $(u,v)$',...
       'fontsize',14,'interpreter','latex','NumColumns',2);







%% Kalman filter from optimal sensors 
S_optimal = S(sort(optimal_input_index,'ascend'),:);

% Measurement noise (Training dataset only)
s_optimal = s(sort(optimal_input_index,'ascend'),:); % optimal sensor inputs 
measurement_noise = s_optimal(:,1:K_training) - S_optimal*a(m,1:K_training);
R_optimal = cov(measurement_noise'); % covariance matrix

% Kalman filter gain L
plant_optimal = ss(A,[B eye(Nx)],S_optimal*C,[],dt_r);
N_optimal = zeros(Nx,length(optimal_input_index)); % noise cross covariance
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

% Plot POD coefficients reconstruction
run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\POD_coeff_plots.m')
% run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\POD_coeff_plots.m');
%clear FIT_a_training FIT_a_validation;
% run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\Reconstruction.m')
% run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\Reconstruction.m')
%save('sysID_kalman_20optimal_sensors','a_e','FIT_a','best_sensors_u','best_sensors_v','best_sensors_uv');