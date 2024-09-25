function [] = initialize_table(filename)

loc= nan;k_in_out_Mad = nan;k_in_out_pMad = nan;isturning_wt = nan;isturning_dKO = nan; ...isturning
Dppc_wt = nan;Dppc_dKO = nan;width_wt = nan;height_wt = nan;width_dKO = nan;height_dKO = nan; ...turning features
phi_pMadn_gsc_wt = nan;phi_pMadcn_precb_wt = nan;phi_pMadn_gsc_dKO = nan;phi_pMadcn_precb_dKO ...phi
 = nan;pMadg_higher = nan;Dad_low = nan;FS_turns_on = nan;higher_branch = nan;FSc_max = nan;perfectness_score = nan; wt_lower_branch = nan;
wt_upper_branch = nan; dadKO_lower_branch = nan; dadKO_upper_branch= nan; ... profile features
pMadg_wt = nan;pMadgn_wt = nan;Dadg_wt = nan;FSg_wt = nan;pMadc_wt = nan;pMadcn_wt = nan;Dadc_wt = nan;FSc_wt = nan; ... ss wt
pMadg_dKO = nan;pMadgn_dKO = nan;pMadc_dKO = nan;pMadcn_dKO = nan;Dadg_dKO = nan;Dadc_dKO = nan;FSg_dKO = nan;FSc_dKO = nan; ... ss dKO
pMadn_LPL_wt = nan;pMadn_LPU_wt = nan;pMadn_LPL_dKO = nan;pMadg_higher = nan; ... bifurcation features
Dpp_LPL_wt = nan;Dpp_LPU_wt = nan;Dpp_LPL_dKO = nan;Dpp_LPU_dKO = nan; ... bf dKO
pMadgn_OScalc = nan;pMadgn_UScalc = nan;pMadgn_OScalctime = nan;pMadgn_UScalctime = nan;... pMadgn response features
pMadgn_OS = nan;pMadgn_US = nan;pMadgn_npeak = nan;pMadgn_Peak = nan;pMadgn_trise = nan;pMadgn_transient = nan;...
pMadgn_setmax = nan;pMadgn_setmin = nan;pMadgn_set = nan;...
pMadcn_OScalc = nan;pMadcn_UScalc = nan;pMadcn_OScalctime = nan;pMadcn_UScalctime = nan;... pMadcn response features
pMadcn_OS = nan;pMadcn_US = nan;pMadcn_npeak = nan;pMadcn_Peak = nan;pMadcn_trise = nan;pMadcn_transient = nan;...
pMadcn_setmax = nan;pMadcn_setmin = nan;pMadcn_set = nan;...
Dadg_OScalc = nan;Dadg_UScalc = nan;Dadg_OScalctime = nan;Dadg_UScalctime = nan;... Dadg RF
Dadg_OS = nan;Dadg_US = nan;Dadg_npeak = nan;Dadg_Peak = nan;Dadg_trise = nan;Dadg_transient = nan;...
Dadg_setmax = nan;Dadg_setmin = nan;Dadg_set = nan;...
Dadc_OScalc = nan;Dadc_UScalc = nan;Dadc_OScalctime = nan;Dadc_UScalctime = nan;... Dadc RF
Dadc_OS = nan;Dadc_US = nan;Dadc_npeak = nan;Dadc_Peak = nan;Dadc_trise = nan;Dadc_transient = nan;...
Dadc_setmax = nan;Dadc_setmin = nan;Dadc_set = nan;...
FSg_OScalc = nan;FSg_UScalc = nan;FSg_OScalctime = nan;FSg_UScalctime = nan;... FSg RF
FSg_OS = nan;FSg_US = nan;FSg_npeak = nan;FSg_Peak = nan;FSg_trise = nan;FSg_transient = nan;...
FSg_setmax = nan;FSg_setmin = nan;FSg_set = nan;...
FSc_OScalc = nan;FSc_UScalc = nan;FSc_OScalctime = nan;FSc_UScalctime = nan;... FSc RF
FSc_OS = nan;FSc_US = nan;FSc_npeak = nan;FSc_Peak = nan;FSc_trise = nan;FSc_transient = nan;...
FSc_setmax = nan;FSc_setmin = nan;FSc_set = nan;...
pMadgn_OScalc_dKO = nan;pMadgn_UScalc_dKO = nan;pMadgn_OScalctime_dKO = nan;pMadgn_UScalctime_dKO = nan;... pMadgn dKO RF
pMadgn_OS_dKO = nan;pMadgn_US_dKO = nan;pMadgn_npeak_dKO = nan;pMadgn_Peak_dKO = nan;pMadgn_trise_dKO = nan;pMadgn_transient_dKO = nan;...
pMadgn_setmax_dKO = nan;pMadgn_setmin_dKO = nan;pMadgn_set_dKO = nan;...
pMadcn_OScalc_dKO = nan;pMadcn_UScalc_dKO = nan;pMadcn_OScalctime_dKO = nan;pMadcn_UScalctime_dKO = nan;... pMadcn dKO RF
pMadcn_OS_dKO = nan;pMadcn_US_dKO = nan;pMadcn_npeak_dKO = nan;pMadcn_Peak_dKO = nan;pMadcn_trise_dKO = nan;pMadcn_transient_dKO = nan;...
pMadcn_setmax_dKO = nan;pMadcn_setmin_dKO = nan;pMadcn_set_dKO = nan;...
Dadg_OScalc_dKO = nan;Dadg_UScalc_dKO = nan;Dadg_OScalctime_dKO = nan;Dadg_UScalctime_dKO = nan;... Dadg dKO RF
Dadg_OS_dKO = nan;Dadg_US_dKO = nan;Dadg_npeak_dKO = nan;Dadg_Peak_dKO = nan;Dadg_trise_dKO = nan;Dadg_transient_dKO = nan;...
Dadg_setmax_dKO = nan;Dadg_setmin_dKO = nan;Dadg_set_dKO = nan;...
Dadc_OScalc_dKO = nan;Dadc_UScalc_dKO = nan;Dadc_OScalctime_dKO = nan;Dadc_UScalctime_dKO = nan;... Dadc dKO RF
Dadc_OS_dKO = nan;Dadc_US_dKO = nan;Dadc_npeak_dKO = nan;Dadc_Peak_dKO = nan;Dadc_trise_dKO = nan;Dadc_transient_dKO = nan;...
Dadc_setmax_dKO = nan;Dadc_setmin_dKO = nan;Dadc_set_dKO = nan;...
FSg_OScalc_dKO = nan;FSg_UScalc_dKO = nan;FSg_OScalctime_dKO = nan;FSg_UScalctime_dKO = nan;... FSg dKO RF
FSg_OS_dKO = nan;FSg_US_dKO = nan;FSg_npeak_dKO = nan;FSg_Peak_dKO = nan;FSg_trise_dKO = nan;FSg_transient_dKO = nan;...
FSg_setmax_dKO = nan;FSg_setmin_dKO = nan;FSg_set_dKO = nan;...
FSc_OScalc_dKO = nan;FSc_UScalc_dKO = nan;FSc_OScalctime_dKO = nan;FSc_UScalctime_dKO = nan;... FSc dKO RF
FSc_OS_dKO = nan;FSc_US_dKO = nan;FSc_npeak_dKO = nan;FSc_Peak_dKO = nan;FSc_trise_dKO = nan;FSc_transient_dKO = nan;...
FSc_setmax_dKO = nan;FSc_setmin_dKO = nan;FSc_set_dKO = nan;...
pMadn_ratio_wt = nan;FS_ratio_wt = nan;pMadn_ratio_dKO = nan;FS_ratio_dKO = nan; datetime = nan; runtime = nan;

writetable(table(loc,k_in_out_Mad,k_in_out_pMad,isturning_wt,isturning_dKO, ...isturning
            Dppc_wt,Dppc_dKO,width_wt,height_wt,width_dKO,height_dKO, ...turning features
            phi_pMadn_gsc_wt,phi_pMadcn_precb_wt,phi_pMadn_gsc_dKO,phi_pMadcn_precb_dKO ...phi
            ,pMadg_higher,Dad_low,FS_turns_on,higher_branch,FSc_max,perfectness_score, ...
            wt_lower_branch,wt_upper_branch,dadKO_lower_branch,dadKO_upper_branch, ... profile features
            pMadg_wt,pMadgn_wt,Dadg_wt,FSg_wt,pMadc_wt,pMadcn_wt,Dadc_wt,FSc_wt, ... ss wt
            pMadg_dKO,pMadgn_dKO,pMadc_dKO,pMadcn_dKO,Dadg_dKO,Dadc_dKO,FSg_dKO,FSc_dKO, ... ss dKO
            pMadn_LPL_wt,pMadn_LPU_wt,pMadn_LPL_dKO,pMadg_higher, ... bifurcation features
            Dpp_LPL_wt,Dpp_LPU_wt,Dpp_LPL_dKO,Dpp_LPU_dKO, ... bf dKO
            pMadgn_OScalc,pMadgn_UScalc,pMadgn_OScalctime,pMadgn_UScalctime,... pMadgn response features
            pMadgn_OS,pMadgn_US,pMadgn_npeak,pMadgn_Peak,pMadgn_trise,pMadgn_transient,...
            pMadgn_setmax,pMadgn_setmin,pMadgn_set,...
            pMadcn_OScalc,pMadcn_UScalc,pMadcn_OScalctime,pMadcn_UScalctime,... pMadcn response features
            pMadcn_OS,pMadcn_US,pMadcn_npeak,pMadcn_Peak,pMadcn_trise,pMadcn_transient,...
            pMadcn_setmax,pMadcn_setmin,pMadcn_set,...
            Dadg_OScalc,Dadg_UScalc,Dadg_OScalctime,Dadg_UScalctime,... Dadg RF
            Dadg_OS,Dadg_US,Dadg_npeak,Dadg_Peak,Dadg_trise,Dadg_transient,...
            Dadg_setmax,Dadg_setmin,Dadg_set,...
            Dadc_OScalc,Dadc_UScalc,Dadc_OScalctime,Dadc_UScalctime,... Dadc RF
            Dadc_OS,Dadc_US,Dadc_npeak,Dadc_Peak,Dadc_trise,Dadc_transient,...
            Dadc_setmax,Dadc_setmin,Dadc_set,...
            FSg_OScalc,FSg_UScalc,FSg_OScalctime,FSg_UScalctime,... FSg RF
            FSg_OS,FSg_US,FSg_npeak,FSg_Peak,FSg_trise,FSg_transient,...
            FSg_setmax,FSg_setmin,FSg_set,...
            FSc_OScalc,FSc_UScalc,FSc_OScalctime,FSc_UScalctime,... FSc RF
            FSc_OS,FSc_US,FSc_npeak,FSc_Peak,FSc_trise,FSc_transient,...
            FSc_setmax,FSc_setmin,FSc_set,...
            pMadgn_OScalc_dKO,pMadgn_UScalc_dKO,pMadgn_OScalctime_dKO,pMadgn_UScalctime_dKO,... pMadgn dKO RF
            pMadgn_OS_dKO,pMadgn_US_dKO,pMadgn_npeak_dKO,pMadgn_Peak_dKO,pMadgn_trise_dKO,pMadgn_transient_dKO,...
            pMadgn_setmax_dKO,pMadgn_setmin_dKO,pMadgn_set_dKO,...
            pMadcn_OScalc_dKO,pMadcn_UScalc_dKO,pMadcn_OScalctime_dKO,pMadcn_UScalctime_dKO,... pMadcn dKO RF
            pMadcn_OS_dKO,pMadcn_US_dKO,pMadcn_npeak_dKO,pMadcn_Peak_dKO,pMadcn_trise_dKO,pMadcn_transient_dKO,...
            pMadcn_setmax_dKO,pMadcn_setmin_dKO,pMadcn_set_dKO,...
            Dadg_OScalc_dKO,Dadg_UScalc_dKO,Dadg_OScalctime_dKO,Dadg_UScalctime_dKO,... Dadg dKO RF
            Dadg_OS_dKO,Dadg_US_dKO,Dadg_npeak_dKO,Dadg_Peak_dKO,Dadg_trise_dKO,Dadg_transient_dKO,...
            Dadg_setmax_dKO,Dadg_setmin_dKO,Dadg_set_dKO,...
            Dadc_OScalc_dKO,Dadc_UScalc_dKO,Dadc_OScalctime_dKO,Dadc_UScalctime_dKO,... Dadc dKO RF
            Dadc_OS_dKO,Dadc_US_dKO,Dadc_npeak_dKO,Dadc_Peak_dKO,Dadc_trise_dKO,Dadc_transient_dKO,...
            Dadc_setmax_dKO,Dadc_setmin_dKO,Dadc_set_dKO,...
            FSg_OScalc_dKO,FSg_UScalc_dKO,FSg_OScalctime_dKO,FSg_UScalctime_dKO,... FSg dKO RF
            FSg_OS_dKO,FSg_US_dKO,FSg_npeak_dKO,FSg_Peak_dKO,FSg_trise_dKO,FSg_transient_dKO,...
            FSg_setmax_dKO,FSg_setmin_dKO,FSg_set_dKO,...
            FSc_OScalc_dKO,FSc_UScalc_dKO,FSc_OScalctime_dKO,FSc_UScalctime_dKO,... FSc dKO RF
            FSc_OS_dKO,FSc_US_dKO,FSc_npeak_dKO,FSc_Peak_dKO,FSc_trise_dKO,FSc_transient_dKO,...
            FSc_setmax_dKO,FSc_setmin_dKO,FSc_set_dKO,...
            pMadn_ratio_wt,FS_ratio_wt,pMadn_ratio_dKO,FS_ratio_dKO,datetime,runtime),filename);

end