% SYS ID & KALMAN FILTER

clear; clc; close all;

%% Load variables

% Vector of spatial locations
load('/home/fs2317/Desktop/FYP/DNS data/x_plotting.mat');
load('/home/fs2317/Desktop/FYP/DNS data/y_plotting.mat');

% Meshgrid
load('/home/fs2317/Desktop/FYP/DNS data/xq.mat');
load('/home/fs2317/Desktop/FYP/DNS data/yq.mat');

% True POD coefficients and eigenmodes
load('/home/fs2317/Desktop/FYP/POD/a.mat');
load('/home/fs2317/Desktop/FYP/POD/phi.mat');

% Mean flow
load('/home/fs2317/Desktop/FYP/POD/Y_mean.mat');

% Computational domain size
x_length = [-2.5 11.5]; y_length = [-4 4];

% Spatial resolution of DNS data
dx = 0.1; dy = 0.1;

% Number of grid points (excluding points inside the cylinders)
n = length(x_plotting);

% Temporal resolution of DNS data
dt = 0.04;
Time = 200:dt:500;
K = length(Time);

%% Sensor library
sensor_library_options = input('Which sensor library do you want to load? \nA == random sampling \nB == POD peaks \nC == single sensors \n');

if sensor_library_options == 'A'
    % ---------------------- Random sampling ----------------------
    step = input('Choose x-y step for the equi-spaced sensor library: ');
    load('/home/fs2317/Desktop/FYP/DNS data/deleteMatrix.mat');
    y1 = 1;
    y2 = 81;
    ym = 0.5*(y2 - y1) + 1;
    y_point = [flip([ym,ym:-step:y1+1,y1]), ym:step:y2-1,y2];
    pos = find(y_point == ym);
    y_point(pos([1,3])) = [];
    
    x1 = 20;
    x2 = 141;
    xm = ceil(0.5*(x2 - x1) + 1);
    x_point = [flip([xm,xm:-step:x1+1,x1]), xm:step:x2-1,x2];
    pos = find(x_point == xm);
    x_point(pos([1,3])) = [];
    
    yqsensor = 0.001*ones(size(yq));
    yqsensor(y_point,[20:step:141]) = yq(y_point,[20:step:141]);
    temp = yqsensor.*deleteMatrix';
    temp(isnan(temp)) = [];
    temp1 = temp(:);
    sensor_point = find(temp1 ~= 0.001);
    % -------------------------------------------------------------
    
elseif sensor_library_options == 'B'
    % ------------------------ POD peaks --------------------------
    number_of_POD_peaks = input('Choose number of POD peaks: ');
    load('/home/fs2317/Desktop/FYP/POD/POD_peaks.mat');
    sensor_point = sort(POD_peaks(1:number_of_POD_peaks),'ascend');
    clear POD_peaks;
    % -------------------------------------------------------------
    
elseif sensor_library_options == 'C'
    % --------------------- Ranked single sensors -----------------
    load('/home/fs2317/Desktop/FYP/Iterative single sensor/simple estimator/best_sensors.mat');
    number_of_best_sensors = input('Choose number of ranked single sensors: ');
    sensor_point = sort(best_sensors(1:number_of_best_sensors),'ascend');
    clear best_sensors;
    % -------------------------------------------------------------
end

plot_sensor_library = input('Do you want to plot the sensor library? (yes=1, no=0): ');
if plot_sensor_library == 1
    figure('name','Sensor library');
    
    if sensor_library_options == 'B'
        plot(x_plotting,y_plotting','k.','Markersize',3); hold on; % mesh points
        plot(x_plotting(sensor_point),y_plotting(sensor_point),'r.','Markersize',18); hold on; % sensors
        for i = 1:length(sensor_point)
            text(x_plotting(sensor_point(i)),y_plotting(sensor_point(i)),[' ',num2str(i)],'Fontsize',12); hold on;
        end
        % cylinders
        pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
        plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
        pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
        plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
        
        grid off; axis equal; view(2);
        xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
        xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
        legend('Mesh','Sensor library','fontsize',18,'interpreter','latex','location','northeastoutside');
        set(gca,'linewidth',1);
        set(gca,'Fontsize',14);
        set(gca,'Ticklength',[0 0]);
        
    else
        plot(x_plotting,y_plotting','k.','Markersize',3); hold on; % mesh points
        plot(x_plotting(sensor_point),y_plotting(sensor_point),'r.','Markersize',18); hold on; % sensors
        
        % cylinders
        pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
        plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
        pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
        plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
        
        grid off; axis equal; view(2);
        xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
        xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
        legend('Mesh','Sensor library','fontsize',18,'interpreter','latex','location','northeastoutside');
        set(gca,'linewidth',1);
        set(gca,'Fontsize',14);
        set(gca,'Ticklength',[0 0]);
    end
end

% Question to proceed/not proceed with current sensor library
sensor_library_proceed = input('Do you want to proceed with this sensor library? (yes=1, no=0): ');
if sensor_library_proceed == 0
    return
end

%% System Identification
m = 1:50; % number of POD modes to estimate
K_training = ceil(K/2); % Training dataset length
K_validation = K - K_training; % Validation dataset length

% n4sid
sysID_operation = input('Have you done sysID already? (yes=1, no=0): ');

if sysID_operation == 0 % Not computed
    
    data = iddata(a(m,:)',[],dt,'Tstart',Time(1)); % create iddata Object for N4SID
    Nx = input('Model order: '); % model order
    opt = n4sidOptions('InitialState','estimate','N4Weight','auto','Focus','Simulation');
    [sysID,x0] = n4sid(data(1:K_training),Nx,'DisturbanceModel','estimate','CovarianceMatrix','estimate',opt); % stochastic identification
    
elseif sysID_operation == 1 % Already computed 
    
    model_order = input('What system order you want to import? ');
    
    switch model_order
        case 'best'
            load('/home/fs2317/Desktop/FYP/Iterative single sensor/kalman estimator/StochasticID_order_best.mat')
        case 50
            load('/home/fs2317/Desktop/FYP/Iterative single sensor/kalman estimator/StochasticID_order_50.mat')
        case 100
            load('/home/fs2317/Desktop/FYP/Iterative single sensor/kalman estimator/StochasticID_order_100.mat')
    end
end

% Plant matrices 
A = sysID.A; % dynamics matrix
B = sysID.B; % input matrix
C = sysID.C; % output matrix   // this is the control law???
Q = sysID.K*sysID.NoiseVariance*sysID.K'; % process noise covariance
Nx = length(x0); % model order

%% Sensor measurements
load('/home/fs2317/Desktop/FYP/POD/Y.mat');
s = zeros(2*length(sensor_point),K);
s(1:2:end,:) = Y(sensor_point,:); % u'
s(2:2:end,:) = Y(sensor_point+n,:); % v'
p = 2*length(sensor_point); % number of p sensor measurements

%% Kalman filter
% Eigenmode topology matrix
S = zeros(p,length(m));
S(1:2:end,:) = phi(sensor_point,m); % u'
S(2:2:end,:) = phi(sensor_point+n,m); % v'

% Measurements noise (Training dataset only)
measurement_noise = s(:,1:K_training) - S*a(m,1:K_training);
R = cov(measurement_noise'); % covariance matrix

% % Cross-covariance (Training dataset only)
% % random process noise vector 
% w = randn(Nx,K_LD); 
% w = chol(Q,'lower')*w; 
% for i = 1:Nx
%     for j = 1:2*length(sensor_point)
%         N(i,j) = (1/K_LD)*w(i,:)*measurement_noise(j,:)';
%     end
% end

% Kalman filter gain L
plant = ss(A,[B eye(Nx)],S*C,[],dt);
N = zeros(Nx,p); % cross-covariance matrix
[kalmf,L,P] = kalman(plant,Q,R,N);

%% Reconstruction of POD coefficients
% Estimation
x_e = x0; % state initialisation
a_e = C*x_e; % pre-allocation
for k = 1:K
    x_e = A*x_e + L*(s(:,k) - S*C*x_e);
    a_e(:,k) = C*x_e;
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

% Plot POD coefficients reconstruction
POD_coefficients_plot;
clear FIT_a_training FIT_a_validation; 

%% Reconstruction of u' and v' fields
velocity_reconstructionFIT;

%% Save results
% name_save = input('Save reconstruction? (if so, type name): ');
% save(name_save,'Nx','p','sensor_point','a_e','FIT_a','FIT_u_training','FIT_u_validation','FIT_v_training','FIT_v_validation','Y_ROM'); 

%% Balanced truncation
kalman_filter = ss(A-L*S*C,L,S*C,[],dt); % kalman filter DT system

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
    dt);

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

figure('name','Optimal sensors');
subplot(5,6,[1:5, 7:11]); % B matrix 
h = imagesc(kalman_filter.B); colormap pink; hold on;
set(gca,'ytick',[]); xticks(sort(optimal_input_index,'ascend')); xtickangle(90);
ad = 0.2*ones(size(kalman_filter.B)); ad(:,optimal_input_index) = 1; set(h,'AlphaData',ad);
xlabel('Number of $2q$ sensor inputs','fontsize',14); ylabel('Number of model states $N_{X}$','fontsize',14); 
title('Optimal $\mathbf{B}$ matrix columns','fontsize',14);

subplot(5,6,[6:6:30]); % s input vector 
h = imagesc(repmat(rand(2*length(sensor_point),1),[1 2])); axis equal; colormap pink; hold on;
set(gca,'xtick',[]); yticks(sort(optimal_input_index,'ascend'));
xlim([0.5 2.5]); ylim([0.5 2*length(sensor_point)+0.5]);
ad = 0.1*ones(2*length(sensor_point),2); ad(optimal_input_index,:) = 1; set(h,'AlphaData',ad);
ylabel('Number of $2q$ sensor inputs','fontsize',14); 
title('Input $\mathbf{s}$','fontsize',14);

%% Best sensors
% Find optimal sensors where | u' | v' | u',v' | are extracted
optimal_sensor_location = ceil(sort(optimal_input_index,'ascend')/2);
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
    
% if any(mod(index_sort,2)) ~= 0
%     optimal_sensor_location_u = optimal_sensor_location_u(optimal_sensor_location_u~=0);
% else
%     optimal_sensor_location_u = [];
% end

% if any(mod(index_sort,2)) ~= 0
%     optimal_sensor_location_v = [];
% else
%     optimal_sensor_location_v = optimal_sensor_location_v(optimal_sensor_location_v~=0);
% end

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
plant_optimal = ss(A,[B eye(Nx)],S_optimal*C,[],dt);
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
POD_coefficients_plot;
clear FIT_a_training FIT_a_validation; 

% save('sysID_kalman_50optimal_sensors','a_e','FIT_a','best_sensors_u','best_sensors_v','best_sensors_uv');