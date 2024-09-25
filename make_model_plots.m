%
clc
clear
close all
%
addpath /scratch/user/razeen/ovary/pss4/Fu_rep/Compute
addpath /scratch/user/razeen/ovary/pss4/Fu_rep/functions

% load('unbiased_screen_parameter_sets.mat')
% load('unbiased_screen_results.mat')
% param_set = sets(:,1:end);

load('unbiased_screen_parameter_sets.mat')
load('directed_screen_results.mat')

param_set = gen_sets(:,1:end);
clear sets

isdadko = 1;
wt = 0;

newwindow = 0;
yesplot = 1;
plotturning = 1;


% loc = 21183; %instance for Fused Repression Model
loc = 2988; % instance for Fused Degrdation Model

sets = param_set(loc,:);

%
% run_palc.m on them
%
data_wt = run_palc(sets);
isturning_wt = data_wt.isturning;

Xwt = data_wt.X;
L_DOT_wt = data_wt.L_DOT;

%
% run_palc.m on dad^KO now for isturning sets
%
data_dKO = run_palc(sets,isdadko);
isturning_dKO = data_dKO.isturning;

XdKO = data_dKO.X;
L_DOT_dKO = data_dKO.L_DOT;

if isturning_wt && isturning_dKO

    %
    % For isturning sets locate_limit_points.m for wt and Dad^KO
    %
    limit_points_wt = locate_limit_points(Xwt);

    Dpp_LPL_wt = limit_points_wt.Dpp_LPL;
    Dpp_LPU_wt = limit_points_wt.Dpp_LPU;
    pMadn_LPL_wt = limit_points_wt.pMad_LPL;
    pMadn_LPU_wt = limit_points_wt.pMad_LPU;

    limit_points_dKO = locate_limit_points(XdKO);
    Dpp_LPL_dKO = limit_points_dKO.Dpp_LPL;
    Dpp_LPU_dKO = limit_points_dKO.Dpp_LPU;
    pMadn_LPL_dKO = limit_points_dKO.pMad_LPL;
    pMadn_LPU_dKO = limit_points_dKO.pMad_LPU;
    %
    % Check the distance between limit points in pMad and Dpp; check if the
    % average of limit points in Dpp is very close to origin
    %
    Dppc_wt = mean([Dpp_LPU_wt,Dpp_LPL_wt]);
    Dppc_dKO = mean([Dpp_LPU_dKO,Dpp_LPL_dKO]);

else
    Dppc_wt = 0.1;
    %leave empty if not turning
    Dpp_LPL_wt = nan;Dpp_LPU_wt = nan;
    pMadn_LPL_wt = nan;pMadn_LPU_wt = nan;
    Dpp_LPL_dKO = nan;Dpp_LPU_dKO = nan;
    pMadn_LPL_dKO = nan;pMadn_LPU_dKO = nan;
    Dppc_dKO = nan;
    width_wt = nan;height_wt = nan;
    width_dKO = nan;height_dKO = nan;
end

%
% At the average of limit points for Dpp simulate_model
%
if yesplot && ~newwindow
    f1 = figure;
    t1 = tiledlayout(1,2);
    ax11 = nexttile;
    axis square
end
sim_model_wt = simulate_model(sets,Dppc_wt,wt,yesplot,newwindow);
pMadcn_wt = sim_model_wt.pMadcn;
FSc_wt = sim_model_wt.FSc;
Tkvactc_wt = sim_model_wt.precb_profile(end,10);

if yesplot && ~newwindow
    ax12 = nexttile;
    axis square
end
sim_model_dKO = simulate_model(sets,Dppc_wt,isdadko,yesplot,newwindow);
pMadcn_dKO = sim_model_dKO.pMadcn;
FSc_dKO = sim_model_dKO.FSc;
Tkvactc_dKO = sim_model_dKO.precb_profile(end,10);

if yesplot
    ax11.YLim(2) = 1.5*max([sim_model_wt.pMadgn sim_model_wt.pMadcn...
        sim_model_dKO.pMadgn sim_model_dKO.pMadcn...
        sim_model_wt.Dadg sim_model_wt.Dadc...
        sim_model_dKO.Dadg sim_model_dKO.Dadc...
        sim_model_wt.FSg sim_model_wt.FSc...
        sim_model_dKO.FSg sim_model_dKO.FSc...
        ]);
    ax11.PlotBoxAspectRatio = [1 1 1];
    ax11.XTick = 0:2:12;
    ax11.YTick = 0:0.2:4;
    ax11.FontName = 'Arial';
    ax11.FontSize = 24;

    ax12.YLim(2) = 1.5*max([sim_model_wt.pMadgn sim_model_wt.pMadcn...
        sim_model_dKO.pMadgn sim_model_dKO.pMadcn...
        sim_model_wt.Dadg sim_model_wt.Dadc...
        sim_model_dKO.Dadg sim_model_dKO.Dadc...
        sim_model_wt.FSg sim_model_wt.FSc...
        sim_model_dKO.FSg sim_model_dKO.FSc...
        ]);
    ax12.PlotBoxAspectRatio = [1 1 1];
    ax12.XTick = 0:2:12;
    ax12.YTick = 0:0.2:4;
    ax12.FontName = 'Arial';
    ax12.FontSize = 24;

    title(t1,['profile# ',num2str(loc)])
    % saveas(gcf,['Plots',filesep,'profile_',num2str(loc),'.png'])
    % close(f1)
end
% exportgraphics(f1,['dynamics',num2str(loc),'.eps'],'ContentType','vector')
% close all

if plotturning && (isturning_dKO && isturning_wt)
    f2 = plot_turning_new(Xwt,XdKO);
    
    
    % % title(['#',num2str(loc)]);
    % f2.Children.XLabel = 'Dpp';
    % f2.Children.YLabel = 'norm. Concentration';
    % f2.Children.XLabel.FontSize = 18;
    % f2.Children.YLabel.FontSize = 18;

    % nexttile(2)
    % hold on
    % plot(Dppc_wt,Tkvactc_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
    % plot(Dppc_wt,Tkvactc_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2) 
    % 
    % nexttile(8)
    % hold on
    % plot(Dppc_wt,FSc_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
    % plot(Dppc_wt,FSc_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2)
    % legend('wt','Dad^{KO}','','','location','best')
    % 
    % nexttile(6)
    % hold on
    % plot(Dppc_wt,pMadcn_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
    % plot(Dppc_wt,pMadcn_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2)
    % exportgraphics(f2,['allbifurcationcurve',num2str(loc),'.jpg'],'ContentType','vector')
    % saveas(gcf,['Plots',filesep,'bifurcation_',num2str(loc),'.png'])
    % close
end
% close all
if plotturning && (isturning_dKO && isturning_wt)
    f3 = figure;
    f3.Position = [680 125 1052 753];
    plot(Xwt(end,:),Xwt(6,:),'k-','linewidth',2.5)
    hold on
    plot(XdKO(end,:),XdKO(6,:),'b-','linewidth',2.5)

    hold on
    plot(Dppc_wt,pMadcn_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
    hold on
    plot(Dppc_wt,pMadcn_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2)

    legend('\it wt','\it dad^{KO}','location','northeast')
    title(['#',num2str(loc)]);
    xlim([0 1])
    ylim([0 3])
    yticks(0:0.2:3)
    xticks(0:0.2:1)
    xlabel('Dpp')
    ylabel('pMadn')
    set(gca,'FontName','Arial','FontSize',24)
    grid on
    grid minor
    axis square
    % exportgraphics(f3,['bifurcationcurve',num2str(loc),'.eps'],'ContentType','vector')
end
% close all
if plotturning && (isturning_dKO && isturning_wt)
    f4 = figure;
    f4.Position = [680 125 1052 753];
    plot(Xwt(end,:),Xwt(6,:),'k-','linewidth',2.5)
    hold on
    plot(Dppc_wt,pMadcn_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)


    legend('\it wt','\it dad^{KO}','location','northeast')
    title(['#',num2str(loc)]);
    xlim([0 1])
    ylim([0 0.25])
    yticks(0:0.05:0.3)
    xticks(0:0.2:1)
    xlabel('Dpp')
    ylabel('pMadn')
    set(gca,'FontName','Arial','FontSize',24)
    grid on
    grid minor
    axis square
    % exportgraphics(f4,['bifurcationcurve_wo_dadko',num2str(loc),'.eps'],'ContentType','vector')
end
%%

if plotturning && (isturning_dKO && isturning_wt)
    f3 = figure;
    f3.Position = [680 125 1052 753];
    plot(Xwt(end,:),Xwt(8,:),'k-','linewidth',2.5)
    hold on
    plot(XdKO(end,:),XdKO(8,:),'b-','linewidth',2.5)

    hold on
    plot(Dppc_wt,FSc_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
    hold on
    plot(Dppc_wt,FSc_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2)

    legend('\it wt','\it dad^{KO}','location','northeast')
    title(['#',num2str(loc)]);
    xlim([0 1])
    ylim([0 1])
    yticks(0:0.2:1)
    xticks(0:0.2:1)
    xlabel('Dpp')
    ylabel('FS')
    set(gca,'FontName','Arial','FontSize',24)
    grid on
    grid minor
    axis square
    % exportgraphics(f3,['FS_bifurcationcurve',num2str(loc),'.eps'],'ContentType','vector')
end
% close all
if plotturning && (isturning_dKO && isturning_wt)
    f4 = figure;
    f4.Position = [680 125 1052 753];
    plot(Xwt(end,:),Xwt(8,:),'k-','linewidth',2.5)
    hold on
    plot(Dppc_wt,FSc_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)


    legend('\it wt','\it dad^{KO}','location','northeast')
    title(['#',num2str(loc)]);
    xlim([0 1])
    ylim([0 1])
    yticks(0:0.2:1)
    xticks(0:0.2:1)
    xlabel('Dpp')
    ylabel('FS')
    set(gca,'FontName','Arial','FontSize',24)
    grid on
    grid minor
    axis square
    % exportgraphics(f4,['FS_bifurcationcurve_wo_dadko',num2str(loc),'.eps'],'ContentType','vector')
end

%%

if plotturning && (isturning_dKO && isturning_wt)
    f3 = figure;
    f3.Position = [680 125 1052 753];
    plot(Xwt(end,:),Xwt(2,:),'k-','linewidth',2.5)
    hold on
    plot(XdKO(end,:),XdKO(2,:),'b-','linewidth',2.5)

    hold on
    plot(Dppc_wt,Tkvactc_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
    hold on
    plot(Dppc_wt,Tkvactc_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2)

    legend('\it wt','\it dad^{KO}','location','northeast')
    title(['#',num2str(loc)]);
    xlim([0 1])
    ylim([0 1])
    yticks(0:0.2:1)
    xticks(0:0.2:1)
    xlabel('Dpp')
    ylabel('Tkv_{act}')
    set(gca,'FontName','Arial','FontSize',24)
    grid on
    grid minor
    axis square
    % exportgraphics(f3,['Tkvact_bifurcationcurve',num2str(loc),'.eps'],'ContentType','vector')
end
% close all
if plotturning && (isturning_dKO && isturning_wt)
    f4 = figure;
    f4.Position = [680 125 1052 753];
    plot(Xwt(end,:),Xwt(2,:),'k-','linewidth',2.5)
    hold on
    plot(Dppc_wt,Tkvactc_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)


    legend('\it wt','\it dad^{KO}','location','northeast')
    title(['#',num2str(loc)]);
    xlim([0 1])
    ylim([0 0.5])
    yticks(0:0.1:0.5)
    xticks(0:0.2:1)
    xlabel('Dpp')
    ylabel('Tkv_{act}')
    set(gca,'FontName','Arial','FontSize',24)
    grid on
    grid minor
    axis square
    % exportgraphics(f4,['Tkvact_bifurcationcurve_wo_dadko',num2str(loc),'.eps'],'ContentType','vector')
end

%%
f5 = figure;
sim_model_wt = simulate_model(sets,Dppc_wt,wt,yesplot,newwindow);

ylim([0 1.5*max([sim_model_wt.pMadgn sim_model_wt.pMadcn...
    sim_model_dKO.pMadgn sim_model_dKO.pMadcn...
    sim_model_wt.Dadg sim_model_wt.Dadc...
    sim_model_dKO.Dadg sim_model_dKO.Dadc...
    sim_model_wt.FSg sim_model_wt.FSc...
    sim_model_dKO.FSg sim_model_dKO.FSc...
    ]) ])
axis square
title(t1,['profile# ',num2str(loc)])
% % exportgraphics(f5,['dynamics_wo_dadko',num2str(loc),'.jpg'],'ContentType','vector')
