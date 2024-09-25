%
clc
clear
close all
%

%
%load parameter sets
%
load("directed_screen_results.mat")
param_set = gen_sets(:,1:end);
clear gen_sets
%
%load results
%
load("directed_screen_results.mat")

%
%Initialize table
%

filename = 'Dppc_screen_Furep.csv';

%
% Grad Sets
%

dpMad = T.pMadgn_wt ./T.pMadcn_wt;
dFS = T.FSc_wt./T.FSg_wt;
dDad = T.Dadg_wt./T.Dadc_wt;

trivial = (dpMad > 1.05) & (dFS > 1.05) & (dDad > 1.05);
clear dpMad dFS dDad
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
clear Dppc_in_preCB Dad_high_Fused_low_GSC pMad_high_Dad_low_CB
tri_bi = trivial & biologically_informed;
T_tri_bi = T(tri_bi,["loc","isturning_wt","Dpp_LPL_wt","Dpp_LPU_wt"]);

%
% Turning Sets
%
bistable = T.isturning_wt == 1;

T_bistable = T(bistable,:);

%
%Grad & BI & turning
%
grad_BI_bistable = trivial & biologically_informed & bistable;
T_grad_BI_bistable = T(grad_BI_bistable,:);


clear T

isdadko = 1;
wt = 0;

newwindow = 0;
yesplot = 0;
plotturning = 1;

for k=1:size(T_grad_BI_bistable,1)

    loc = T_grad_BI_bistable.loc(k);
    sets = param_set(loc,:);

    %
    % access palc data for the k'th set
    %
    isturning_wt = T_grad_BI_bistable.isturning_wt(k);
    isturning_dKO = T_grad_BI_bistable.isturning_dKO(k);



    if isturning_wt && isturning_dKO
        %
        % For isturning sets access limit points for wt and Dad^KO
        %

        Dpp_LPL_wt = T_grad_BI_bistable.Dpp_LPL_wt(k);
        Dpp_LPU_wt = T_grad_BI_bistable.Dpp_LPU_wt(k);
        pMadn_LPL_wt = T_grad_BI_bistable.pMadn_LPL_wt(k);
        pMadn_LPU_wt = T_grad_BI_bistable.pMadn_LPU_wt(k);
        pMadn_LPL_dKO = T_grad_BI_bistable.pMadn_LPL_dKO(k);

        data_dKO = run_palc(sets,isdadko);
        
        XdKO = data_dKO.X;
        limit_points_dKO = locate_limit_points(XdKO);
        pMadn_LPU_dKO = limit_points_dKO.pMad_LPU;
        
        %
        % Check the distance between limit points in pMad and Dpp; check if the
        % average of limit points in Dpp is very close to origin
        %
        Dppc_wt = mean([Dpp_LPU_wt,Dpp_LPL_wt]);

    else

        continue;

    end
    
    Dpp_left = linspace(Dpp_LPU_wt,Dppc_wt,5);
    Dpp_right = linspace(Dppc_wt,Dpp_LPL_wt,5);
    Dppc = [Dpp_left Dppc_wt Dpp_right];

   if ~issorted(Dppc)
    continue
   end


    parfor i = 1:length(Dppc)

        sim_model_wt(i) = simulate_model(sets,Dppc(i),wt,yesplot,newwindow);
        % FSC_wt(k,i) = sim_model_wt.FSc;

        sim_model_dKO(i) = simulate_model(sets,Dppc(i),isdadko,yesplot,newwindow);


    end

data(k).sim_model_wt = sim_model_wt;
data(k).sim_model_dKO = sim_model_dKO;
disp(["Done = ",k])
end
save('pDppcScreen.mat','data','-v7.3')
