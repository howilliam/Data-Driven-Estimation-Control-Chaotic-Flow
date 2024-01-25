% SCRIPT TO EXTRACT DATA FROM CFD SIMULATIONS 

clear; clc; close all;

%% Mesh
x_length = [-3 11]; y_length = [-4 4]; % size of flow domain
dx = 0.1; dy = 0.1; % discretisation
x_grid = x_length(1):dx:x_length(2); 
y_grid = y_length(1):dy:y_length(2);
[xq,yq] = meshgrid(x_grid,y_grid); 

%% Exclude points inside the geometry
excludex = [ones(1,(-0.5-x_length(1))/dx), zeros(1,(0.5+0.5)/dx+1), ones(1,round((x_length(2)-0.5)/dx))];
excludey = [ones(1,(-1.5-y_length(1))/dy), zeros(1,round((-0.5+1.5)/dy+1)), ones(1,round((0.5+0.5)/dy-1)), zeros(1,round((1.5-0.5)/dy+1)), ones(1,(y_length(2)-1.5)/dy)]';

deleteMatrix = ones(length(excludex), length(excludey)); % pre-allocation
for i = 1:length(excludex)
    for j = 1:length(excludey)
        if excludex(i) == 0 && excludey(j) == 0
            deleteMatrix(i,j) = NaN;
        end
    end
end
clear excludex excludey;

%% Assemble snapshot matrices 
flow_time = 20.01:0.01:100;

% read only mesh coordinates
name_file = ['XYZ_Internal_Table_table_',num2str(flow_time(2)*100,'%3.2f'),'.1'];
fid = fopen(name_file);
fgetl(fid); % remove 1st line of text
get_data = textscan(fid,'%*f %f %f %*f %*f'); 
x = get_data{1,1}; y = get_data{1,2}; 
fclose(fid);

% pre-allocation of velocity field
u = zeros(11259,1); v = zeros(11259,1); 
interpolation_u = griddata(x,y,u,xq,yq).*deleteMatrix'; interpolation_u(isnan(interpolation_u)) = [];
interpolation_v = griddata(x,y,v,xq,yq).*deleteMatrix'; interpolation_v(isnan(interpolation_v)) = [];
U(:,length(flow_time)) = interpolation_u(:); 
V(:,length(flow_time)) = interpolation_v(:); 
% end pre-allocation

processes = 4; % number of parallel processes
parfor(n = 1:length(flow_time))
    
    name_file = ['FFF-3-',num2str(flow_time(n),'%3.2f'),'.1'];
    fid = fopen(name_file);
    fgetl(fid); % remove 1st line of text
        if n == 1
        %                        node number x-coord y-coord u     v  x-coord y-coord
        get_data = textscan(fid,'%*f         %*f     %*f     %f    %f %*f     %*f');
        u = get_data{1,1}; v = get_data{1,2};
        fclose(fid);   
        else
        %                        node number x-coord y-coord u     v 
        get_data = textscan(fid,'%*f         %*f     %*f     %f    %f');
        u = get_data{1,1}; v = get_data{1,2};
        fclose(fid);
        end 
    
        interpolation_u = griddata(x,y,u,xq,yq).*deleteMatrix'; interpolation_u(isnan(interpolation_u)) = [];
        interpolation_v = griddata(x,y,v,xq,yq).*deleteMatrix'; interpolation_v(isnan(interpolation_v)) = [];
        U(:,n) = interpolation_u(:);
        V(:,n) = interpolation_v(:);
end
interpolation_x = xq.*deleteMatrix'; interpolation_x(isnan(interpolation_x)) = [];
interpolation_y = yq.*deleteMatrix'; interpolation_y(isnan(interpolation_y)) = []; 
x_plotting = interpolation_x(:); 
y_plotting = interpolation_y(:); 

%% Save variables of interest 
save('U.mat','U'); save('V.mat','V'); 
save('x_plotting','x_plotting'); save('y_plotting','y_plotting'); 
save('xq','xq'); save('yq','yq'); 
save('deleteMatrix','deleteMatrix');