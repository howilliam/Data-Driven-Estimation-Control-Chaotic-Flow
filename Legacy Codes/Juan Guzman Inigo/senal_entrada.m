%%%%% system identification input %%%%%

%% input signal %%

N=5000;                                      %Nombre de pas de temps
fny=1/(2*5);                                %Frequence de Nyquist 1/(2*dt)
fphys= 0.012;                              %Frequence du ph�nomene physique
finput=20*fphys/fny;                   %frequence de l'input normalis�e, 20 x plus grand que la physique
u=idinput(N,'rbs',[0 finput], [-1 1]);    %vecteur des input en fction du temps

plot(u);

% sampling=1/fny/2;
% N=length(u);
% n=0:N/2;
% 
% fn=n/N/sampling;
% 
% fid = fopen('rbs.txt','w');
% fprintf(fid,'%12i\n',length(u));
% for k=1:length(u)
% fprintf(fid,'%12.8f\n',u(k));
% end
% fclose(fid);
% figure,
% yk=fft(u);
% semilogy(fn,2/N*abs(yk(1:length(fn))),'-'),
% %plot(fn,2/N*abs(yk(1:length(fn))),'-'),

shg,
