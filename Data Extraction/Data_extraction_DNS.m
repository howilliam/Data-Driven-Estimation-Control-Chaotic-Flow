%% DNS DATA EXTRACTION 
% This is the first codes from a series for the 2022 Dissertation project
% by Juan Antonio Auli Gomez
% This code extracts the data from the flow domain stored in csv files, an
% output from star-ccm+'s xyz internal tables. 
clear; clc;


%path='\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\HPC data\DNS DATA (mat
%files)'
% path='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\ExtractedData\';  %HPC path
% pathSAVE='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\Save Test\'; %the path it saves to 
% pathDNS='C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\mat files\'; %generates mat files ???
%path='/home/wh219/FYP/ExtractedData/'  %HPC path open loop
path='/home/wh219/FYP/Controller/Data Extraction Controller and Excitation/'  %HPC path exciatation / controller
%pathSAVE='/home/wh219/FYP/Open Loop/' %path to save the grid data open loop
pathSAVE='/home/wh219/FYP/Controller/Controller Save/' %path to save the grid data excitation / controller
pathDNS='/home/wh219/FYP/DNS mat files/'; % This is the file path to save the generated mat files ???

%% Presentation grid definition
    % This is the domain used in STAR-CCM+ where the probes were placed
    
x_length=[-3 11];
y_length=[-4 4]; % Domain coordinates for sensor grid in star-ccm
% the nondimensional spacing for each is 0.1, defined as x/l and y/h
dx= 0.1; dy=0.1; %% from length/ no: probes == 14/140 & 8/80

x_grid= [x_length(1):dx:x_length(2)]; % creates the set of points in the grispace. These match the points in the imported csv files 
y_grid= [y_length(1):dy:y_length(2)];

% create a meshgrid for 2d plotting 
[xq,yq]= meshgrid(x_grid,y_grid);

%% File name definition for reading data

ft=10010:10:100000;   % flow time // file time
%creates vector with time intervals at which data was saved in the DNS 
% % digits(6);
% % flow_time=vpa(ft,6); % obtains number of decimal places in ft such that it matches the CSVfiles
filetime = arrayfun(@string, ft); %creates a string from the syms above

%% SNAPSHOT MATRIX ASSEMBLY

prompt= input('Do you want to import GRID data from CSV files?')
if prompt==1
    %[U V P]=get_data(path,filetime);
    [U, V]=get_data(path,filetime);
    %% Sort sensors by location
           
        % Import x y locations of grid sensors 
        [X_sensors]=readmatrix([path 'XYZ_Internal_Table_table_20050.csv'],'Range','C:C'); % x coordinates of simulation
        [Y_sensors]=readmatrix([path 'XYZ_Internal_Table_table_20050.csv'], 'Range','D:D'); % y coordinates of simulation


        %spatial_coords
        % x_sort=[X_sensors;X_inGeo];
        % y_sort= [Y_sensors;Y_inGeo];

        % sort the matrices
            Sort_U2=[X_sensors Y_sensors U]; 
            Sort_V2=[X_sensors Y_sensors V];
            %Sort_P=[X_sensors Y_sensors P];
            Sort_U=sortrows(Sort_U2,[1,2]);
            Sort_V=sortrows(Sort_V2,[1,2]);
            %Sort_P=sortrows(Sort_P,[1,2]);
        % extract the organised data

            x_plotting=Sort_U(1:end,1);
            y_plotting=Sort_U(1:end,2);
            
            
            U=Sort_U(1:end,3:end);
            V=Sort_V(1:end,3:end);
            %P=Sort_P(:,3:end);
            %clear Sort_P Sort_V Sort_U;
            clear Sort_V Sort_U;
            
            %ASK TO SAVE THE DATA 
            if input('would you like to save the GRID data?')==1

            %     name= prompt('enter filename');
              % save([pathSAVE 'Griddata_HPC_200500_002TR20.mat'],'U','V','P','x_plotting','y_plotting','xq','yq');
               save([pathSAVE 'Griddata_Control_9000_03_05.mat'],'U','V','x_plotting','y_plotting','xq','yq');
               %save([pathSAVE 'Testing.mat'],'U','V','x_plotting','y_plotting','xq','yq');
            end
            
            
else 
    prompt=input('Do you want to import GRID RAW data fom saved . mat files?')
     if prompt==1 
         
         load([path 'UVP_RAW_24jun.mat']);

         
         
     elseif input('Do you want to import GRID SORTED data from saved in .mat files?')==1

        %Import sorted data
        
         load([path 'Griddata_HPC_002TR40.mat']);
       
     end
end



%% Obtain fish data 

% if input('Import fish data from CSV files?')==1
%     fishdata
% elseif input('Import data from .mat files')==1
%          load([path 'fish_sorted_HPC_200500_002TR20.mat']);
% 
% end 