
%% QR pivoting 

%Plant is loaded before
%Eigenmodes, will use S= Phi_f*H

%number of POD MODES
p_modes=length(m);

%% QR pivoting

    % DT Lyapunov equation for the controllability Gramian 
    Wc = dlyap(sysID.A,sysID.K*(sysID.K)');
    
    % Output controllability Gramian
    Woc = H*sysID.C*Wc*(H*sysID.C)';
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
    [~,R_qr,pivot] = qr((Phi_fish(:,1:p_modes)*Lc)','vector'); 
    
    % Save optimally ranked sensors
     %QRsensors = pivot(1:p_modes)'; %will have same number of sensors as modes
    QRsensors = pivot(1:10)'; %will have same number of sensors as modes
    
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
    S_optimal = phi(sort(QRsensors,'ascend'),m)*H(m,:);
    % Measurement noise (Training dataset only)
%     s_optimal = s(sort(QRsensors,'ascend'),:); % optimal sensor inputs 
    s_optimal = Y_mat_fish(sort(QRsensors,'ascend'),:); % optimal sensor inputs 
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
index_v2 = QRsensors_v - 11259;   
% Plot POD coefficients reconstruction
run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\FIT_U_V.m')
run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\Reconstructionv2.m')
run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\POD_coeff_plots.m')
% run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\POD_coeff_plots.m');
%clear FIT_a_training FIT_a_validation;
%run('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\MSCproj\200500 data\supporting\Reconstruction.m')
% run('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\supporting\Reconstruction.m')
%save('sysID_kalman_20optimal_sensors','a_e','FIT_a','best_sensors_u','best_sensors_v','best_sensors_uv');
    