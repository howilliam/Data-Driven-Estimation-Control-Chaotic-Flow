%% SVD fish sensors
clear;clc;
% function [a_P,lam_p,Phi]=SVD_fish(P_fish,)
%%Temporal definitions
dt=0.002; %temporal resolution of data extraction
r=40; % rate at which the data is extracted
T=100020*dt:dt*r:250220*dt;% flow time
K=length(T); %number of snapshots

%% LOAD DATASETS
pathSVD= '/home/jaa21/Downloads/Project/26 jul test/';

DATASET = 'A'


if DATASET=='A'
    load([pathSVD 'fish_sorted_HPC_200500_002TR20.mat'])
elseif DATASET=='B'
    load([pathSVD 'fish_sortedB_HPC_200500_002TR20.mat'])
else 
    error('Error: no dataset selected')
end 

%% Perform SVD for fish data
    
    MODES=5;
    snapshot=1;
    %Pressure data
   [a_P,Y_P,Y_P_mean,lam_p,Phi_P]=SVDfish(P_fish,snapshot,MODES,X_fish,Y_fish,'Pressure');
    
    % Shear Stress data
    [a_WS,Y_WS,Y_WS_mean,lam_WS,Phi_WS]=SVDfish(WS_fish,snapshot,MODES,X_fish,Y_fish, 'Shear Stress');
    
    % Save the pod data
    %pathSVD='/home/jaa21/Downloads/Project/POD files/'
    %save([pathSVD 'fishsvd_200500_002TR40.mat'], 'a_P','a_WS','Y_P','Y_WS','Y_P_mean','Y_WS_mean','lam_p','lam_WS','Phi_P','Phi_WS');
  
