% QR PIVOTING OPTIMIZATION VALIDATION 

clear; clc; close all; 
set(0,'defaulttextinterpreter', 'latex');
set(0,'DefaultTextFontname', 'latex');
set(0,'DefaultAxesFontName', 'latex');
set(0,'defaultAxesFontSize',14);

% Load the optimisation history 
load('/home/fs2317/Desktop/FYP/Optimal sensor placement/modes 1-11 (27 sensors and POD input)/optimisationBT.mat');

% Extract information 
N = 2*length(optimisationBT.sensor_point); % number of sensor inputs 
class = 2; % set the number of sensors to test in combination 
comb = factorial(N)/(factorial(class)*factorial(N - class)); % number of combinations 
c = subsetcombs(1:N, comb, class); % call function

% Find effectiveness of optimisation 
B = optimisationBT.model_original.B;
Wo_b = optimisationBT.Wo_b; 

obj = zeros(1,length(comb)); obj2 = obj; 
for i = 1:comb
    B_star = B(:,c(i,:)); % create sensor-group actuator matrix
    obj(i) = logdet(B_star'*Wo_b*B_star); % log determinant  
end 

% Optimal sensors
S_B = optimisationBT.S_B; 
S_B = S_B(:,1:2); 
B_star = B*S_B;
obj_opt = logdet(B_star'*Wo_b*B_star);

effectiveness = length(find(obj_opt > obj))/comb*100; % [%]

% Plot
figure('name','Effectiveness of QR optimisation');
subplot(1,3,[1 2]);
sorted = sort(obj,'ascend');
b = bar([sorted, obj_opt],1,'FaceColor','flat'); hold on; 
x = [comb-40, comb+1, comb+1, comb-40, comb-40]; y = [sorted(comb-40), sorted(comb-40), max(sorted), max(sorted), sorted(comb-40)];
plot(x,y,'k-','LineWidth',2);
b.CData(end,:) = [1 0 0];
xlabel('Number of combinations','fontsize',14,'interpreter','latex');
ylim([min(obj) max(obj)]); xlim([1 comb]);
title(['$\log \left | \mathbf{B^{*^{\top}} \widehat{\mathbf{W}}_0 \: \mathbf{B^{*}} } \right |$ for all combinations of ',num2str(N),...
       ' sensor  inputs grouped in classes of ',num2str(class)],'fontsize',16,'interpreter','latex');

subplot(1,3,3);
b = bar([sorted(end-30:end), obj_opt],1,'FaceColor','flat'); hold on; 
b.CData(end,:) = [1 0 0];
yl = yline(obj_opt,'r--',['> ',num2str(round(effectiveness,2)),'% of all combinations'],'LineWidth',2);
yl.LabelHorizontalAlignment = 'left'; yl.FontSize = 14;
ylim([min(sorted(end-30:end)) max(sorted(end-30:end))]); xticks([]);


function c = subsetcombs(v, k, p)
   %v vector to pick from
   %k number of combinations to return
   %p number of elements to pick in each combination
     N = numel(v);
     assert(k <= nchoosek(N, p), 'You''ve requested more unique combinations than exist');
     c = zeros(k, p);
     for cidx = 1:k
        while true
           c(cidx, :) = v(randperm(N, p));  %generate random combination
           if ~any(arrayfun(@(row) isempty(setdiff(c(row, :), c(cidx, :))), 1:cidx-1))
              %no duplicate, move onto next cidx
              break;
           end
           %duplicate combination, try another one
        end
     end
end
  
function v = logdet(A, op)
%LOGDET Computation of logarithm of determinant of a matrix
%
%   v = logdet(A);
%       computes the logarithm of determinant of A. 
%
%       Here, A should be a square matrix of double or single class.
%       If A is singular, it will returns -inf.
%
%       Theoretically, this function should be functionally 
%       equivalent to log(det(A)). However, it avoids the 
%       overflow/underflow problems that are likely to 
%       happen when applying det to large matrices.
%
%       The key idea is based on the mathematical fact that
%       the determinant of a triangular matrix equals the
%       product of its diagonal elements. Hence, the matrix's
%       log-determinant is equal to the sum of their logarithm
%       values. By keeping all computations in log-scale, the
%       problem of underflow/overflow caused by product of 
%       many numbers can be effectively circumvented.
%
%       The implementation is based on LU factorization.
%
%   v = logdet(A, 'chol');
%       If A is positive definite, you can tell the function 
%       to use Cholesky factorization to accomplish the task 
%       using this syntax, which is substantially more efficient
%       for positive definite matrix. 
%
%   Remarks
%   -------
%       logarithm of determinant of a matrix widely occurs in the 
%       context of multivariate statistics. The log-pdf, entropy, 
%       and divergence of Gaussian distribution typically comprises 
%       a term in form of log-determinant. This function might be 
%       useful there, especially in a high-dimensional space.       
%
%       Theoretially, LU, QR can both do the job. However, LU 
%       factorization is substantially faster. So, for generic
%       matrix, LU factorization is adopted. 
%
%       For positive definite matrices, such as covariance matrices,
%       Cholesky factorization is typically more efficient. And it
%       is STRONGLY RECOMMENDED that you use the chol (2nd syntax above) 
%       when you are sure that you are dealing with a positive definite
%       matrix.
%
%   Examples
%   --------
%       % compute the log-determinant of a generic matrix
%       A = rand(1000);
%       v = logdet(A);
%
%       % compute the log-determinant of a positive-definite matrix
%       A = rand(1000);
%       C = A * A';     % this makes C positive definite
%       v = logdet(C, 'chol');
%

%   Copyright 2008, Dahua Lin, MIT
%   Email: dhlin@mit.edu
%
%   This file can be freely modified or distributed for any kind of 
%   purposes.
%

% argument checking

assert(isfloat(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
    'logdet:invalidarg', ...
    'A should be a square matrix of double or single class.');

if nargin < 2
    use_chol = 0;
else
    assert(strcmpi(op, 'chol'), ...
        'logdet:invalidarg', ...
        'The second argument can only be a string ''chol'' if it is specified.');
    use_chol = 1;
end

% computation

if use_chol
    v = 2 * sum(log(diag(chol(A))));
else
    [L, U, P] = lu(A);
    du = diag(U);
    c = det(P) * prod(sign(du));
    v = log(c) + sum(log(abs(du)));
end
end 