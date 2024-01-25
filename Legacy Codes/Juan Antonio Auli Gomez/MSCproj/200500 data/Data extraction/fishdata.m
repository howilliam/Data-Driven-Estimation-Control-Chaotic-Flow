%% fish data collection 
%this code extracts the data from the fish sensors specified in star-ccm+
%S
%ft=7040:40:25000;
%% Import data 
    %import xyz coordinates from sensors
    [X_fish]=readmatrix([path 'FISHsensors_t_7540.csv'],'Range','F:F');
    [Y_fish]=readmatrix([path 'FISHsensors_t_7540.csv'], 'Range','G:G');
    [Z_fish]=readmatrix([path 'FISHsensors_t_7540.csv'],'Range','H:H');
    XYZ_RAW=[X_fish Y_fish Z_fish];
    % import Wall shear stress and presure data
    processes=15;
    parfor n=1:length(ft)
        name_file=[path 'FISHsensors_t_', cell2mat(filetime(n)), '.csv'];
        Pfish(:,n)= readmatrix(name_file, 'Range','A:A');
        WS(:,n)= readmatrix(name_file, 'Range','B:B');
    end


%% sort the data 

    % Create a matrix with xzy values for sorting
    xyz=[Z_fish X_fish Y_fish];
    Pxy = [xyz Pfish];
    WSxy= [xyz WS];
    %discard negative z values and remove z values from dataset (2:end)
    Pxy=Pxy(find(Pxy(:,1)>0),2:end);
    WSxy=WSxy(find(WSxy(:,1)>0),2:end);
    %Remove the doubble layer data 
%     PLEASE READ THE README, THE BELOW SECTION IS ONLY APPICABLE FOR THE
%     CURRENT PROJECT, THE FOLLOWING PROCESS ONLY ENSURES THAT THE DATA
%     CORRESPONDING TO THE WALL VLUES ARE KEPT. HOWEVER THIS CHANGES BASED
%     ON YOUR GEOMETRY. 

%   BE ADVISED THAT IT IS BETTER TO USE THE CODES AS A GUIDE AND ONLY TO
%   ADA

    Psort=[];
    WSsort=[];
    i=[30 68 165 314 375 384 467 559 595 657 766 919]; e=[59 75 254 371 376 388 543 573 616 695 836 1001];
    
    for g=1:length(i)
        Psort=[Psort; Pxy(i(g):e(g),:)];
        WSsort=[WSsort; WSxy(i(g):e(g),:)];
               
    end; %clear g xyz Pxy WSxy %Pfish WS;

    
    % Sort by Z values (twice created by STAR-CCm+)
%     WSsort=sortrows(WSsort,[1 2]); 
%     Psort=sortrows(Psort,[1 2]);
    
    %% %save the value of missing sensors
    
    % Organize the data by upper and lower sorted by x
    upperWS=WSsort(find(WSsort(:,2)>10.001),:);
    upperWS=sortrows(upperWS,1);
    lowerWS=WSsort(find(WSsort(:,2)<=10.001),:);
    lowerWS=sortrows(lowerWS,1);
    
    upperP=Psort(find(Psort(:,2)>10.001),:);
    upperP=sortrows(upperP,1);
    lowerP=WSsort(find(Psort(:,2)<=10.001),:);
    lowerP=sortrows(lowerP,1);
    
    
    %% REASSEMBLE
%     WSsort=[upperWS; lowerWS]; %<<<<<<<<<<<<unocomment to sort by U-L
%     Psort=[upperP; lowerP];
    WSsort=zeros(size(WSsort));
    Psort=zeros(size(Psort));
    
    WSsort(1:2:end,:)=upperWS;
    WSsort(2:2:end,:)=lowerWS;
    
    
    Psort(1:2:end,:)=upperP;
    Psort(2:2:end,:)=lowerP;
    %save xy for future plotting
    %f_u=upperP(:,1:2);f_l=lowerP(:,1:2);
    %save([path 'F_U_L_GEO.mat']);
    %clear upperP lowerP upperWS lowerWS
    
    %Split the data from x y coordinates
    X_fish = WSsort(:,1);
    Y_fish = WSsort(:,2);
    WS_fish=WSsort(:,3:end);
    P_fish=Psort(:,3:end);
    %clear Psort WSsort
%% save the sorted data
if input('save the sorted FISH data?')==1    
    save([pathSAVE 'fish_sorted_HPC_200500_002TR20.mat'],'P_fish','WS_fish','X_fish','Y_fish');
end

