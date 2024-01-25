%%


load('\\icnas1.cc.ic.ac.uk\jaa21\Desktop\MSCproj\200500 data\fish_sorted_HPC_200500_002TR20.mat','Y_fish','X_fish');% contains data order as used in codes
plot(X_fish([1:2:500 500]),Y_fish([1:2:500 500]),'k', 'linewidth', 1.3); hold on;
plot(X_fish([1 2:2:500]),Y_fish([1 2:2:500]),'k', 'linewidth', 1.3); axis equal
ylim([9.5 10.5]); 
xlim([7.5 12.5])
title(['Fish, Top 20 ' plotname ' sensors (by average)'],'interpreter','latex','fontsize',15)

s=scatter(X_fish(TOP_20(:,1)),Y_fish(TOP_20(:,1)),'ro')
%labelpoints(X_fish(TOP_20(:,1)),Y_fish(TOP_20(:,1)),[1:20],'NE')

% insert datatip
%         
%         dt = datatip(s,'DataIndex',10,'interpreter','latex','fontsize',10);
%         s.DataTipTemplate.DataTipRows(1).Label = '$St$'; s.DataTipTemplate.DataTipRows(1).Format = '%0.2f';
%         s.DataTipTemplate.DataTipRows(2).Label = '$|a(St)|$'; s.DataTipTemplate.DataTipRows(2).Format = '%0.2f';