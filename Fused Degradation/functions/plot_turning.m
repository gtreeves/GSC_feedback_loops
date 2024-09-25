function f = plot_turning(Xwt,XdadKO)
%PLOT_TURNING Summary of this function goes here
%   Detailed explanation goes here
f = figure;
f.WindowState = 'maximized';
t = tiledlayout(2,4);
titlestring = {'Tkv','Tkv^{act}','Mad','Madn','pMad','pMadn','Dad','FS'};

for j = 1:size(Xwt,1)-1
    nexttile
    plot(Xwt(end,:),Xwt(j,:),'k-','linewidth',2.5)
    hold on
    plot(XdadKO(end,:),XdadKO(j,:),'b-','linewidth',2.5)
%         hold on
%         plot(XrfxKO(10,:),XrfxKO(j,:),'r-','linewidth',1.5)
%         hold on
%         plot(Xrfxi(10,:),Xrfxi(j,:),'g--','linewidth',1.5)
    hold off
    xticks(0:0.2:1)
    title(titlestring(j))
    set(gca,'FontSize',14)
end
legend('wt','Dad^{KO}','location','best')
set(gca,'FontSize',14)
end

