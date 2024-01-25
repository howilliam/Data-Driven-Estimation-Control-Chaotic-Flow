function [POD_peaks] = maxPOD_sensor(phi,x_plotting,m)
%maxPOD_sensor: calculates the point of max POD spatial mode for
%               identification of optimal sensor location

% INPUTS: 
% phi = POD spatial modes matrix
% x_plotting = x coordinate vector
% m = POD modes [1,2,3,...m]

% OUTPUTS:
% POD_peaks = spatial coordinate of max spatial POD mode

%% Max spatial POD
n = length(x_plotting); % number of spatial locations

for i = 1:length(m) 
    PHI = (phi(1:n,m(i)).^2 + phi(n+1:end,m(i)).^2).^(0.5);
    POD_peaks(i,:) = find(PHI == max(PHI));
end

end