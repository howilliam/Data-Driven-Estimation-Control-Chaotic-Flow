%RECONSTRUCTION VIDEO/PLOT CODE 
pathSID = '/home/wh219/FYP/Controller/Controller Save/';
m = 1:50;
K = 9000;
load('/home/wh219/FYP/Open Loop/POD_Open_9000.mat');
Y_open = Y;
Y_mean_open = Y_mean;
phi_open = phi;
load('/home/wh219/FYP/Open Loop/9000_SYSID_Nx_50_M_50.mat');
x0_open = x0;
A_open = sysID.A; % dynamics matrix
C_open = sysID.C; %output matrix

sensor_point = [15442 15420 16171 16230 4327 6449 20620 3674 22157 18987 7086 17617 3778 10490 4334 7803 21967 18203 3300 5567];
p = length(sensor_point);
s = zeros(p,K);
s = Y_open(sensor_point,:); 

S = zeros(p,length(m));
S = phi_open(sensor_point,m); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

x_e = x0_open; % state initialisation
a_e = C_open*x_e; % pre-allocation
for k = 1:K
    x_e = A_open*x_e + L*(s(:,k) - S*C_open*x_e); % Remember this S includes Phi_fish*H
    a_e(:,k) = C_open*x_e;
end
load('/home/wh219/FYP/Controller/Controller Save/Griddata_Excited.mat','y_plotting','x_plotting');% contains data order as used in codes
load('/home/wh219/FYP/Controller/Controller Save/POD_Control_v2_9000.mat');
load([pathSID 'Excite_9000_v2SYSID_Nx50_M_50.mat'] );
x0_CTRL = x0;
A_CTRL = sysID.A; % dynamics matrix
B_CTRL = sysID.B; %input matrix
C_CTRL = sysID.C; %output matrix

s2 = zeros(p,K);
s2 = Y(sensor_point,:); 

S2 = zeros(p,length(m));
S2 = phi(sensor_point,m); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
     m_interest = 1:10; % MODES OF INTEREST
     %m_interest=m;
%% define and load necesary data for plot
%load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\cyl_geo.mat')
%load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\fishpoints_ORIG.mat')

x_length=[-3 11];
y_length=[-4 4]; % Domain coordinates for sensor grid in star-ccm
dx= 0.1; dy=0.1;
n = length(x_plotting); % number of grid points (excluding points inside the cylinders)
DNSvsREC = figure('Name','DNS vs. Reconstruction: velocity magnitude');


%snapshot = find(T == 15000*0.004); % choose snapshot to compare
%snapshot=1200;



folder = ['Controller Recording']
pathvid=['/home/wh219/FYP/Controller/Controller Recording/']; % FOLDER WHERE THE VIDEO WILL BE SAVED
mkdir(pathvid) % IF FOLDER DOES NOT EXIST, IT WILL CREATE THE FOLDER


            
            
            obj=VideoWriter([pathvid 'Control_v2.avi'],'Motion JPEG AVI');
                
   

            obj.FrameRate=10;
            obj.Quality=100;

            open(obj);
            fr=1;
            Y_REC = zeros(size(Y_mean)); % pre-allocation
            Y_REC = Y_REC + phi_open(:,1:50)*a_e(1:50,:);
            for snapshot=100:200  %% juan had 3756 diferent time
            %for snapshot=2700:2800  %% juan had 3756 diferent time
            Magnitude_DNS = sqrt((Y_open(1:n,snapshot) + Y_mean_open(1:n)).^2 + (Y_open(n+1:end,snapshot) + Y_mean_open(n+1:end)).^2); % DNS
            Magnitude_CTRL = sqrt((Y(1:n,snapshot) + Y_mean(1:n)).^2 + (Y(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % DNS Control
            
            %SNAPSHOT TO COMPARE


            %% plot
          
%                 DNS_contours = subplot(3,5,[1 2]);
                DNS_contours = subplot(2,1,1);
                zi = griddata(x_plotting,y_plotting,Magnitude_DNS,xq,yq);
                contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on; % the 4th index number of contours, the more contours, the smoother
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); 
                colormap(DNS_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 2.5]);
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title(['DNS (timestep=' num2str(snapshot) ')'],'fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
                
                CTRL_contours = subplot(2,1,2);
                zi = griddata(x_plotting,y_plotting,Magnitude_CTRL,xq,yq);
                contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); 
                colormap(CTRL_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 2.5]);
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title(['Control (timestep=' num2str(snapshot) ')'],'fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
               
     
                F(fr)=getframe(DNSvsREC);
                fr=fr+1;


            end
            writeVideo(obj,F)
            obj.close();




%%

