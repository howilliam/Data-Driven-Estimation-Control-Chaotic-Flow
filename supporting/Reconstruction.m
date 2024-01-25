%RECONSTRUCTION VIDEO/PLOT CODE 


%load('C:\Users\devil\OneDrive - Imperial College London\Documents\FYP\Save Test\Griddata_HPC_200500_002TR20.mat','y_plotting','x_plotting');% contains data order as used in codes
%   load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\400500 data\Griddata_HPC_400500_002TR20.mat');  
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



folder = ['Nx ' num2str(Nx) 'm ' num2str(m(end))]
pathvid=['/home/wh219/FYP' folder '/']; % FOLDER WHERE THE VIDEO WILL BE SAVED
mkdir(pathvid) % IF FOLDER DOES NOT EXIST, IT WILL CREATE THE FOLDER

if input('Create Video?')==1

            switch type
                case 1
                    obj=VideoWriter([pathvid 'NX_' num2str(Nx) 'm_' num2str(m(end)) 'sens_' num2str(length(sensor_point)) '_' plotname ],'MPEG-4');
                case 2
                    obj=VideoWriter([pathvid 'NX_' num2str(Nx) 'm_' num2str(m(end)) 'sens_' num2str(r) '_optimal sesnors'],'MPEG-4');
            end 

            obj.FrameRate=10;
            obj.Quality=100;

            open(obj);
            fr=1;
            for snapshot=1500:1600  %% juan had 3756 diferent time
            %for snapshot=2700:2800  %% juan had 3756 diferent time
            Magnitude_DNS = sqrt((Y(1:n,snapshot) + Y_mean(1:n)).^2 + (Y(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % DNS

            %SNAPSHOT TO COMPARE


            %% plot
                Y_RECONSTRUCTED = zeros(size(Y_mean)); % pre-allocation
                for i = 1:length(m_interest)
                    Y_RECONSTRUCTED = Y_RECONSTRUCTED + phi(:,m_interest(i))*a_e(m_interest(i),:); 
                end
                Magnitude_ROM = sqrt((Y_RECONSTRUCTED(1:n,snapshot) + Y_mean(1:n)).^2 + (Y_RECONSTRUCTED(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % ROM
                
%                 DNS_contours = subplot(3,5,[1 2]);
                DNS_contours = subplot(2,1,1);
                zi = griddata(x_plotting,y_plotting,Magnitude_DNS,xq,yq);
                contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); 
                colormap(DNS_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title(['DNS (timestep=' num2str(snapshot) ')'],'fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
                %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');

%                 ESTIMATED_contours = subplot(3,5,[4 5]); 
                ESTIMATED_contours = subplot(2,1,2); 
                zi = griddata(x_plotting,y_plotting,Magnitude_ROM,xq,yq);
                contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); 
                colormap(ESTIMATED_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
                s=scatter(x_plotting(index_u2),y_plotting(index_u2),'filled','rs');  %% u velocity index being plotted for sensors
                %s=scatter(x_plotting(index(6:end)),y_plotting(index(6:end)),'bo');
                hold on;
    
                s=scatter(x_plotting(index_v2),y_plotting(index_v2),'filled','bo');    %% v velocity index being plotted for sensors
                hold on
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title(['Reconstruction'],'fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
                %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');

%                 FITTING = subplot(3,5,[6:15]);
%                 FIT = 100*(1 - abs(Magnitude_ROM - Magnitude_DNS)./abs(Magnitude_DNS));
%                 zi = griddata(x_plotting,y_plotting,FIT,xq,yq);
%                 zi = max(zi,0); % remove all negative fit values
%                 contourf(xq,yq,zi,80,'Edgecolor','flat'); hold on;
%                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
%                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
%                 pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
%                 plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
%                 shading interp; axis equal; view(2); colormap(FITTING,jet); c = colorbar('eastoutside'); 
%                 caxis([0 100]); colorbar('Ticks',[0:25:100]);
%                 xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
%                 xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
%                 title('FIT [\%]','fontsize',14,'interpreter','latex');
%                 set(gca,'linewidth',1);
%                 set(gca,'TickLength',[0 0]);
%                 set(gca,'Fontsize',12); 
                
                
                
                F(fr)=getframe(DNSvsREC);
                fr=fr+1;


            end
            writeVideo(obj,F)
            obj.close();


elseif input('create snapshot in time?')
                
            snapshot=300;
            Magnitude_DNS = sqrt((Y(1:n,snapshot) + Y_mean(1:n)).^2 + (Y(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % DNS

            %SNAPSHOT TO COMPARE


            % plot
                Y_RECONSTRUCTED = zeros(size(Y_mean)); % pre-allocation
                for i = 1:length(m_interest)
                    Y_RECONSTRUCTED = Y_RECONSTRUCTED + phi(:,m_interest(i))*a_e(m_interest(i),:); 
                end
                Magnitude_ROM = sqrt((Y_RECONSTRUCTED(1:n,snapshot) + Y_mean(1:n)).^2 + (Y_RECONSTRUCTED(n+1:end,snapshot) + Y_mean(n+1:end)).^2); % ROM

                DNS_contours = subplot(3,5,[1 2]);
                zi = griddata(x_plotting,y_plotting,Magnitude_DNS,xq,yq);
                contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); 
                colormap(DNS_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title(['DNS (timestep=' num2str(snapshot) ')'],'fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
                %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');

                ESTIMATED_contours = subplot(3,5,[4 5]); 
                zi = griddata(x_plotting,y_plotting,Magnitude_ROM,xq,yq);
                contourf(xq,yq,zi,50,'Edgecolor','flat'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); 
                colormap(ESTIMATED_contours,othercolor('RdBu10')); colorbar('eastoutside'); caxis([0 max(Magnitude_DNS)]);
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title(['Reconstruction'],'fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
                %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');

                FITTING = subplot(3,5,[6:15]);
                FIT = 100*(1 - abs(Magnitude_ROM - Magnitude_DNS)./abs(Magnitude_DNS));
                zi = griddata(x_plotting,y_plotting,FIT,xq,yq);
                zi = max(zi,0); % remove all negative fit values
                contourf(xq,yq,zi,80,'Edgecolor','flat'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[1.5 0.5 0.5 1.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none'); hold on;
                pgon = polyshape([-0.5 -0.5 0.5 0.5],[-0.5 -1.5 -1.5 -0.5]);
                plot(pgon,'FaceColor',[.5 .5 .5],'FaceAlpha',1,'EdgeColor','none');
                shading interp; axis equal; view(2); colormap(FITTING,jet); c = colorbar('eastoutside'); 
                caxis([0 100]); colorbar('Ticks',[0:25:100]);
                xlim([xq(1,1) xq(1,end)]); ylim([yq(1,1) yq(end,1)]);
                xticks(0:2:round(x_length(2))); yticks(y_length(1):2:y_length(2));
                title('FIT [\%]','fontsize',14,'interpreter','latex');
                set(gca,'linewidth',1);
                set(gca,'TickLength',[0 0]);
                set(gca,'Fontsize',12); 
                %fill(s_geo,h_geo_u,'k'); hold on; fill(s_geo,h_geo_L,'k');hold on;fill(x_cyl,y_cyl,'k');
                
                
                saveas(gcf,[pathvid 'ReconstructionModes_' num2str(j) '_' num2str(length(sensor_point))],'epsc')

else 
    
end

%%



