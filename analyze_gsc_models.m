clc
clear
close all 

%% Load the correct version
%For either Fused repression or degradation model
%options: 'unbiased_screen_results' OR 'directed_screen_results'

% load('unbiased_screen_parameter_sets.mat')
% load('unbiased_screen_results.mat')

load('unbiased_screen_parameter_sets.mat')
load('directed_screen_results.mat')
sets = gen_sets;

%% Categories

%
% Trivial Sets
%

dpMad = T.pMadgn_wt ./T.pMadcn_wt;
dFS = T.FSc_wt./T.FSg_wt;
dDad = T.Dadg_wt./T.Dadc_wt;

trivial = (dpMad > 1.05) & (dFS > 1.05) & (dDad > 1.05);
T_trivial = T(trivial,:);
%% what are trivial sets

T1 = T(~trivial,:);
trivialLocs = T1.loc;
setsTrivialLoc = sets(trivialLocs,:);

%
% Biologically-informed parameter sets
%

% location of Dppc in turning sets (what about width/height)
Dppc_in_preCB = T.Dppc_wt >= 0.05;
% In the GSC, concentration of Dad > Fused
Dad_high_Fused_low_GSC = T.Dadg_wt > T.FSg_wt;
% In the CB, pMad > Dad
pMad_high_Dad_low_CB = T.pMadc_wt > T.Dadc_wt;
% Dad should fall to nearly zero in CB (Dad_low = Dadc_wt < 0.05)
Dad_low = T.Dad_low;
% Fused levels in CB should be higher than GSC (FS_turns_on = FSc_wt > 1.5*FSg_wt;)
Fused_turns_on = T.FS_turns_on;

biologically_informed = Dppc_in_preCB & Dad_high_Fused_low_GSC & pMad_high_Dad_low_CB & Dad_low & Fused_turns_on;
T_biologically_informed = T(biologically_informed,:);
%
% Bistable Sets
%
bistable = T.isturning_wt == 1;

T_bistable = T(bistable,:);

%
%Grad and BI
%
trivial_BI = trivial & biologically_informed;
T_trivial_BI = T(trivial_BI,:);

%
%Grad & BI & turning
%
trivial_BI_bistable = trivial & biologically_informed & bistable;
T_trivial_BI_bistable = T(trivial_BI_bistable,:);

%
% Grad turning
%
trivial_bistable = bistable & trivial;
T_trivial_bistable = T(trivial_bistable,:);
%% Figures for manuscript

% 1. Distribution of u11, u12, mu_Tkv and KFS

parameter_list = {'u_{1}','u_{2}','u_{3}','u_{4}','u_{5}','u_{6}','u_{7}','u_{8}','u_{9}','u_{10}','u_{11}','u_{12}','\mu_{Tkv}','\mu_{Mad}','K_{Dad}','K_{FS}','\tau_{Tkv}','\tau_{Mad}','\tau_{Dad}','\tau_{FS}'};

plot_sets_grad = sets(T_trivial.loc,:);
plot_sets_biologically_informed = sets(T_biologically_informed.loc,:);
plot_sets_turning = sets(T_bistable.loc,:);
plot_sets_grad_bi = sets(T_trivial_BI.loc,:);
plot_sets_grad_bi_turning = sets(T_trivial_BI_bistable.loc,:);
% tiledlayout('flow')
for i = [11,12,13,16]
    F1 =figure;
    
    F1.Units = 'inches';
    F1.Position(3:4) = [1.5 1.5]; %W,H
    histogram(log10(plot_sets_grad_bi(:,i)),'Normalization','probability','BinWidth',0.25,'FaceColor',"#4DBEEE",'FaceAlpha',0.5)
    hold on
    histogram(log10(plot_sets_grad_bi_turning(:,i)),'Normalization','probability','BinWidth',0.25,'FaceColor',"#EDB120",'FaceAlpha',0.5)
    hold off
    if i == 16
        xticks([-2 -1 0 1])
        xlim([-2 1])
    else
        xticks([-2 -1 0 1 2])
        xlim([-2 2])
    end
    axis square
    set(gca,'FontSize',10,'FontName','Arial')
    exportgraphics(F1,['figure',num2str(i),'.eps'],'ContentType','vector')
    close all
end
%% sensitivity coefficient and delta-sensitivity coefficient

A = T_trivial_BI.phi_pMadn_gsc_wt;
B = T_trivial_BI.phi_pMadn_gsc_dKO;
C = sum((B-A) > 0)/size(T_trivial_BI,1);
close all
f31 = figure;
f31.Position = [680 125 1052 753];
% t3 = tiledlayout(1,2);
% nexttile
clear c
% figure
c(1) = cdfplot(T_trivial_BI.phi_pMadn_gsc_wt);
hold on
c(2) = cdfplot(T_trivial_BI.phi_pMadn_gsc_dKO);
xlim([0 1])
title('')
lgd = legend('\it wt','\it dad^{KO}','location','best');
xlabel('\phi')
ylabel('CDF')
xticks(0:0.25:1)
axis square
set(c,'LineWidth',2,{'Color'},{'k','k'}',{'LineStyle'},{'-','--'}')
set(gca,'FontSize',24,'FontName','Arial')
% exportgraphics(f31,'sensitivity.eps','ContentType','vector')
% delphi grad bi sets
close all
clear c
f32 = figure;
f32.Position = [680 125 1052 753];
% nexttile
c(1) = cdfplot(T_trivial_BI.phi_pMadn_gsc_wt - T_trivial_BI.phi_pMadn_gsc_dKO);
xlim([-0.5 1])
% title(lgd,['— GSC',' -- preCB'])
xlabel('\phi_{wt} - \phi_{dad^{KO}}')
ylabel('CDF')
title('')
xticks(-1:0.25:1)
axis square
set(c,'LineWidth',2,{'Color'},{'k'}',{'LineStyle'},{'-'}')
set(gca,'FontSize',24,'FontName','Arial')

% exportgraphics(f32,'dsensitivity.eps','ContentType','vector')

%% response speed

close all
clear c
f4 = figure;
f4.Position = [680 125 1052 753];
% tiledlayout(1,2)
% nexttile
c(1) = cdfplot(T_trivial_BI.FSc_trise/60);
hold on
c(2) = cdfplot(T_trivial_BI.FSc_trise_dKO/60);

lgd = legend('\it wt','\it dad^{KO}','location','southeast');
title(['Fu/Smurf response time for Biologically informed sets'])
% title('Dad ratio')
xlabel('time [min]')
ylabel('CDF')
% xlim([0 1])
xticks(0:100:800)
axis('square')
set(c,'LineWidth',2,{'Color'},{'k','k'}',{'LineStyle'},{'-','--'}')
set(gca,'FontSize',24,'FontName','Arial')
exportgraphics(f4,'responsetime.eps','ContentType','vector')

% %'s

N_grad = round(100 * size(T_trivial,1)/size(T,1));

N_grad_bi = round(100 * size(T_trivial_BI,1)/size(T,1));

N_grad_bi_turning = round(100 * size(T_trivial_BI_bistable,1)/size(T_trivial_BI,1));

N_Dad_low = round(100 * sum(T_bistable.Dad_low == 1)/(size(T_bistable.wt_lower_branch,1)));

N_FS_turns_on = round(100 * sum(T_bistable.Dad_low == 1)/(size(T_bistable.wt_lower_branch,1)));

N_pMad_lower_branch_wt = round(100 * sum(T_trivial_BI_bistable.wt_lower_branch == 1)/sum(T_trivial_BI_bistable.wt_lower_branch == 1 | T_trivial_BI_bistable.wt_upper_branch == 1));
E = round(100 * sum(T_trivial_BI_bistable.wt_lower_branch == 1)/(sum(T_trivial_BI_bistable.wt_lower_branch == 1) + sum(T_trivial_BI_bistable.wt_upper_branch == 1)));
F = round(100 * sum(T_trivial_BI_bistable.wt_upper_branch == 1)/size(T_trivial_BI_bistable,1));
N_pMad_lower_branch_wt_ = round(100 * sum(T_bistable.wt_lower_branch == 1)/sum(T_bistable.wt_lower_branch == 1 | T_bistable.wt_upper_branch == 1));

N_pMad_dadko_higher_branch = round(100 * sum(T_trivial_BI_bistable.higher_branch == 1)/sum(T_trivial_BI_bistable.wt_lower_branch == 1));
%% piechart grad-bi/all

figure
pie([size(T_trivial_BI,1) size(T,1)-size(T_trivial_BI,1)])
exportgraphics(gca,'gradbi_furep.eps','ContentType','vector','Append',true)
exportgraphics(gca,[model,'.pdf'],'ContentType','vector','Append',true)


%% Histograms of parameter distributions

% load('screen_sets_Fu_deg.mat')

parameter_list = {'u_{1}','u_{2}','u_{3}','u_{4}','u_{5}','u_{6}','u_{7}','u_{8}','u_{9}','u_{10}','u_{11}','u_{12}','u_{13}','\mu_{Tkv}','\mu_{Mad}','K_{Dad}','K_{FS}','\tau_{Tkv}','\tau_{Mad}','\tau_{Dad}','\tau_{FS}'};

plot_sets_grad = sets(T_trivial.loc,:);
plot_sets_grad_bi = sets(T_trivial_BI.loc,:);
plot_sets_turning = sets(T_bistable.loc,:);

f1 = figure;
% f1 = figure;
f1.Units = 'inches';
f1.Position(3:4) = [6 5]; %W,H
% f1.WindowState = "fullscreen";
tiledlayout('flow')
for i = 1:21
    nexttile
    % histogram(log10(plot_sets_grad(:,i)),'Normalization','probability','NumBins',10,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.5)
    hold on
    histogram(log10(plot_sets_grad_bi(:,i)),'Normalization','probability','NumBins',10,'FaceColor',"#4DBEEE",'FaceAlpha',0.5)
    histogram(log10(plot_sets_turning(:,i)),'Normalization','probability','NumBins',10,'FaceColor',"#EDB120",'FaceAlpha',0.5)
    hold off
    title(parameter_list(i))
    set(gca,'FontSize',8,'FontName','Arial')
end
legend('Biologically-informed','Bistable')

% exportgraphics(f1,[model,'.pdf'],'ContentType','vector','Append',true)
exportgraphics(f1,'allParameters.eps','ContentType','vector')