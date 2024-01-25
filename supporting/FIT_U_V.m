% Reconstruction of u' and v' fields
Y_ROM = zeros(size(Y_mean)); % pre-allocation
n = length(x_plotting);
for i = 1:length(m)
    Y_ROM = Y_ROM + phi(:,m(i))*a_e(m(i),:);
end

% FIT [%]
FIT_u_training = 100*(1-((norm(Y(1:n,1:K_training) - Y_ROM(1:n,1:K_training)))/(norm(Y(1:n,1:K_training) - mean(Y(1:n,1:K_training))))));
FIT_u_validation = 100*(1-((norm(Y(1:n,K_training+1:end) - Y_ROM(1:n,K_training+1:end)))/(norm(Y(1:n,K_training+1:end) - mean(Y(1:n,K_training+1:end))))));
FIT_v_training = 100*(1-((norm(Y(n+1:end,1:K_training) - Y_ROM(n+1:end,1:K_training)))/(norm(Y(n+1:end,1:K_training) - mean(Y(n+1:end,1:K_training))))));
FIT_v_validation = 100*(1-((norm(Y(n+1:end,K_training+1:end) - Y_ROM(n+1:end,K_training+1:end)))/(norm(Y(n+1:end,K_training+1:end) - mean(Y(n+1:end,K_training+1:end))))));