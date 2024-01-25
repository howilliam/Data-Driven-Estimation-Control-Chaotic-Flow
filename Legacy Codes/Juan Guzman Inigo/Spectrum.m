%======================================
%
%    Frequency spectrum of a signal
%
%======================================

% Step 1 of code to find fphys for the senal_entrada code

close all,
file='POD/sensor.dat';
NSENSOR = 4;

 % Select sensor position.
               % 1 x=100, 2 x=150, 3 x=200
               % 4 x=250, 5 x=300, 6 x=350
               % 7 x=400, 8 x=450, 9 x=500
               % 10 x=600, 11 x=700, 12 x=750
               % 13 x=800, 14 x=850

data=load(file); % so data are my pod coefficients at sensor points 1-15 and 4000 takes in all the time
[N,M]=size(data); % N = 4000 , M = 15   %% for me i need 1800 x 50 or 10 so in my a (my coefficients)
t=data(:,1);  % 4000
m=data(:,NSENSOR); % m = 4000 x 1   takes the fourth sensor ???

%plot(t,m)
[Y0,f0]=fourier_transform(t,m);
% m=data(:,3);
% [Y1,f1]=fourier_transform(t,m);
% m=data(:,5);
% [Y2,f2]=fourier_transform(t,m);
% m=data(:,8);
% [Y3,f3]=fourier_transform(t,m);

semilogy(f0,Y0,'c'),xlabel('Frequency(Hz)'),ylabel('|Y|')
% hold on,
% semilogx(f1,Y1,'g'),xlabel('Frequency(Hz)'),ylabel('|Y|')
% hold on,
% semilogx(f2,Y2,'r'),xlabel('Frequency(Hz)'),ylabel('|Y|')
% hold on,
% semilogx(f3,Y3),xlabel('Frequency(Hz)'),ylabel('|Y|'),legend('2','3','5','8')
% hold off,

%file='w.dat';
%data=load(file);
%w=data(:,2)/0.25;
%var(w)

%alpha=0.9;
%h=kstest(w)
%=======================================