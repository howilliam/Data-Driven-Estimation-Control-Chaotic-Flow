%% DATA EXTRACTION
%function [U V P]=get_data(path,filetime)
function [U, V]=get_data(path,filetime)
processes = 15; %number of parallel processes
%parfor (n= 1:length(ft))
%for n= 645:length(ft) 
parfor n=1:length(filetime)
    name_file=[path 'XYZ_Internal_Table_table_', cell2mat(filetime(n)), '.csv'];
    
    U(:,n)= readmatrix(name_file, 'Range','A:A');
    V(:,n)= readmatrix(name_file, 'Range','B:B');
    %P(:,n)= readmatrix(name_file, 'Range','E:E'); %UNCOMENT IF PRESSURE IS
    %NEEDED
       
end



