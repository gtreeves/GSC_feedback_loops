%
clc
clear
close all


%Unbiased screen
load("unbiased_screen_parameter_sets.mat")
param_set = sets(:,1:end);
clear sets

%Directed screen
% load("directed_screen_parameter_sets.mat")
% param_set = gen_sets(:,1:end);
% clear gen_sets

total_sets = size(param_set,1);

isdadko = 1;
wt = 0;

newwindow = 0;
yesplot = 1;
plotturning = 1;
delta = 1;

if yesplot && ~exist('Plots','dir')
    mkdir('Plots')
end

filename = ['Fu_rep',datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.csv'];
initialize_table(filename);

parfor i=1:total_sets
     tic

        sets = param_set(i,:);
        %

        %
        % is kin > kout
        %
        k_in_out_Mad = sets(:,8)/sets(:,7);
        k_in_out_pMad = sets(:,4)/sets(:,5);
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
    
            width_wt = Dpp_LPL_wt - Dpp_LPU_wt;
            height_wt = pMadn_LPU_wt - pMadn_LPL_wt;
    
            width_dKO = Dpp_LPL_dKO - Dpp_LPU_dKO;
            height_dKO = pMadn_LPU_dKO - pMadn_LPL_dKO;
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
        end
        sim_model_wt = simulate_model(sets,Dppc_wt,wt,yesplot,newwindow);

        pMadg_wt = sim_model_wt.pMadg;
        pMadgn_wt = sim_model_wt.pMadgn;
        Dadg_wt = sim_model_wt.Dadg;
        FSg_wt = sim_model_wt.FSg;        

        pMadc_wt = sim_model_wt.pMadc;
        pMadcn_wt = sim_model_wt.pMadcn;
        Dadc_wt = sim_model_wt.Dadc;
        FSc_wt = sim_model_wt.FSc;

        pMadn_ratio_wt = pMadgn_wt/pMadcn_wt;
        FS_ratio_wt = FSc_wt/FSg_wt;
% Get response features from simulate model wt for pMadn, Dad and FS
%         gsc_wt = sim_model_wt.gsc_profile;
%         precb_wt = sim_model_wt.precb_profile;
        pMadgn_OScalc = sim_model_wt.response_features(6).OS_calc;
        pMadgn_UScalc = sim_model_wt.response_features(6).US_calc;
        pMadgn_OScalctime = sim_model_wt.response_features(6).OS_calc_time;
        pMadgn_UScalctime = sim_model_wt.response_features(6).US_calc_time;  
        pMadgn_OS = sim_model_wt.response_features(6).Overshoot;
        pMadgn_US = sim_model_wt.response_features(6).Undershoot;
        pMadgn_npeak = length(sim_model_wt.response_features(6).pks);
        pMadgn_Peak = sim_model_wt.response_features(6).Peak;
        pMadgn_trise = sim_model_wt.response_features(6).RiseTime;
        pMadgn_transient = sim_model_wt.response_features(6).TransientTime;
        pMadgn_setmax = sim_model_wt.response_features(6).SettlingMax;
        pMadgn_setmin = sim_model_wt.response_features(6).SettlingMin;
        pMadgn_set = sim_model_wt.response_features(6).SettlingTime;

        Dadg_OScalc = sim_model_wt.response_features(7).OS_calc;
        Dadg_UScalc = sim_model_wt.response_features(7).US_calc;
        Dadg_OScalctime = sim_model_wt.response_features(7).OS_calc_time;
        Dadg_UScalctime = sim_model_wt.response_features(7).US_calc_time;        
        Dadg_OS = sim_model_wt.response_features(7).Overshoot;
        Dadg_US = sim_model_wt.response_features(7).Undershoot;
        Dadg_npeak = length(sim_model_wt.response_features(7).pks);
        Dadg_Peak = sim_model_wt.response_features(7).Peak;
        Dadg_trise = sim_model_wt.response_features(7).RiseTime;
        Dadg_transient = sim_model_wt.response_features(7).TransientTime;
        Dadg_setmax = sim_model_wt.response_features(7).SettlingMax;
        Dadg_setmin = sim_model_wt.response_features(7).SettlingMin;
        Dadg_set = sim_model_wt.response_features(7).SettlingTime;
        Dadg_peaktime = sim_model_wt.response_features(7).PeakTime;


        FSg_OScalc = sim_model_wt.response_features(8).OS_calc;
        FSg_UScalc = sim_model_wt.response_features(8).US_calc;
        FSg_OScalctime = sim_model_wt.response_features(8).OS_calc_time;
        FSg_UScalctime = sim_model_wt.response_features(8).US_calc_time;        
        FSg_OS = sim_model_wt.response_features(8).Overshoot;
        FSg_US = sim_model_wt.response_features(8).Undershoot;
        FSg_npeak = length(sim_model_wt.response_features(8).pks);
        FSg_Peak = sim_model_wt.response_features(8).Peak;
        FSg_trise = sim_model_wt.response_features(8).RiseTime;
        FSg_transient = sim_model_wt.response_features(8).TransientTime;
        FSg_setmax = sim_model_wt.response_features(8).SettlingMax;
        FSg_setmin = sim_model_wt.response_features(8).SettlingMin;
        FSg_set = sim_model_wt.response_features(8).SettlingTime;
        FSg_peaktime = sim_model_wt.response_features(8).PeakTime;

        pMadcn_OScalc = sim_model_wt.response_features(14).OS_calc;
        pMadcn_UScalc = sim_model_wt.response_features(14).US_calc;
        pMadcn_OScalctime = sim_model_wt.response_features(14).OS_calc_time;
        pMadcn_UScalctime = sim_model_wt.response_features(14).US_calc_time;  
        pMadcn_OS = sim_model_wt.response_features(14).Overshoot;
        pMadcn_US = sim_model_wt.response_features(14).Undershoot;
        pMadcn_npeak = length(sim_model_wt.response_features(14).pks);
        pMadcn_Peak = sim_model_wt.response_features(14).Peak;
        pMadcn_trise = sim_model_wt.response_features(14).RiseTime;
        pMadcn_transient = sim_model_wt.response_features(14).TransientTime;
        pMadcn_setmax = sim_model_wt.response_features(14).SettlingMax;
        pMadcn_setmin = sim_model_wt.response_features(14).SettlingMin;
        pMadcn_set = sim_model_wt.response_features(14).SettlingTime;

        Dadc_OScalc = sim_model_wt.response_features(15).OS_calc;
        Dadc_UScalc = sim_model_wt.response_features(15).US_calc;
        Dadc_OScalctime = sim_model_wt.response_features(15).OS_calc_time;
        Dadc_UScalctime = sim_model_wt.response_features(15).US_calc_time; 
        Dadc_OS = sim_model_wt.response_features(15).Overshoot;
        Dadc_US = sim_model_wt.response_features(15).Undershoot;
        Dadc_npeak = length(sim_model_wt.response_features(15).pks);
        Dadc_Peak = sim_model_wt.response_features(15).Peak;
        Dadc_trise = sim_model_wt.response_features(15).RiseTime;
        Dadc_transient = sim_model_wt.response_features(15).TransientTime;
        Dadc_setmax = sim_model_wt.response_features(15).SettlingMax;
        Dadc_setmin = sim_model_wt.response_features(15).SettlingMin;
        Dadc_set = sim_model_wt.response_features(15).SettlingTime;
        Dadc_peaktime = sim_model_wt.response_features(15).PeakTime; 

        FSc_OScalc = sim_model_wt.response_features(16).OS_calc;
        FSc_UScalc = sim_model_wt.response_features(16).US_calc;
        FSc_OScalctime = sim_model_wt.response_features(16).OS_calc_time;
        FSc_UScalctime = sim_model_wt.response_features(16).US_calc_time; 
        FSc_OS = sim_model_wt.response_features(16).Overshoot;
        FSc_US = sim_model_wt.response_features(16).Undershoot;
        FSc_npeak = length(sim_model_wt.response_features(16).pks);
        FSc_Peak = sim_model_wt.response_features(16).Peak;
        FSc_trise = sim_model_wt.response_features(16).RiseTime;
        FSc_transient = sim_model_wt.response_features(16).TransientTime;
        FSc_setmax = sim_model_wt.response_features(16).SettlingMax;
        FSc_setmin = sim_model_wt.response_features(16).SettlingMin;
        FSc_set = sim_model_wt.response_features(16).SettlingTime;
        FSc_peaktime = sim_model_wt.response_features(16).PeakTime;        
        %
        % At the average of limit points for Dpp simulate_model for dad^KO
        %
        if yesplot && ~newwindow
            ax12 = nexttile;
        end
        sim_model_dKO = simulate_model(sets,Dppc_wt,isdadko,yesplot,newwindow);
        pMadg_dKO = sim_model_dKO.pMadg;
        pMadgn_dKO = sim_model_dKO.pMadgn;
        pMadc_dKO = sim_model_dKO.pMadc;
        pMadcn_dKO = sim_model_dKO.pMadcn;

        Dadg_dKO = sim_model_dKO.Dadg;
        Dadc_dKO = sim_model_dKO.Dadc;
        
        FSg_dKO = sim_model_dKO.FSg;
        FSc_dKO = sim_model_dKO.FSc;

        pMadn_ratio_dKO = pMadgn_dKO/pMadcn_dKO;
        FS_ratio_dKO = FSc_dKO/FSg_dKO;

        pMadgn_OScalc_dKO = sim_model_dKO.response_features(6).OS_calc;
        pMadgn_UScalc_dKO = sim_model_dKO.response_features(6).US_calc;
        pMadgn_OScalctime_dKO = sim_model_dKO.response_features(6).OS_calc_time;
        pMadgn_UScalctime_dKO = sim_model_dKO.response_features(6).US_calc_time;  
        pMadgn_OS_dKO = sim_model_dKO.response_features(6).Overshoot;
        pMadgn_US_dKO = sim_model_dKO.response_features(6).Undershoot;
        pMadgn_npeak_dKO = length(sim_model_dKO.response_features(6).pks);
        pMadgn_Peak_dKO = sim_model_dKO.response_features(6).Peak;
        pMadgn_trise_dKO = sim_model_dKO.response_features(6).RiseTime;
        pMadgn_transient_dKO = sim_model_dKO.response_features(6).TransientTime;
        pMadgn_setmax_dKO = sim_model_dKO.response_features(6).SettlingMax;
        pMadgn_setmin_dKO = sim_model_dKO.response_features(6).SettlingMin;
        pMadgn_set_dKO = sim_model_dKO.response_features(6).SettlingTime;

        Dadg_OScalc_dKO = sim_model_dKO.response_features(7).OS_calc;
        Dadg_UScalc_dKO = sim_model_dKO.response_features(7).US_calc;
        Dadg_OScalctime_dKO = sim_model_dKO.response_features(7).OS_calc_time;
        Dadg_UScalctime_dKO = sim_model_dKO.response_features(7).US_calc_time;        
        Dadg_OS_dKO = sim_model_dKO.response_features(7).Overshoot;
        Dadg_US_dKO = sim_model_dKO.response_features(7).Undershoot;
        Dadg_npeak_dKO = length(sim_model_dKO.response_features(7).pks);
        Dadg_Peak_dKO = sim_model_dKO.response_features(7).Peak;
        Dadg_trise_dKO = sim_model_dKO.response_features(7).RiseTime;
        Dadg_transient_dKO = sim_model_dKO.response_features(7).TransientTime;
        Dadg_setmax_dKO = sim_model_dKO.response_features(7).SettlingMax;
        Dadg_setmin_dKO = sim_model_dKO.response_features(7).SettlingMin;
        Dadg_set_dKO = sim_model_dKO.response_features(7).SettlingTime;
        Dadg_peaktime_dKO = sim_model_dKO.response_features(7).PeakTime;


        FSg_OScalc_dKO = sim_model_dKO.response_features(8).OS_calc;
        FSg_UScalc_dKO = sim_model_dKO.response_features(8).US_calc;
        FSg_OScalctime_dKO = sim_model_dKO.response_features(8).OS_calc_time;
        FSg_UScalctime_dKO = sim_model_dKO.response_features(8).US_calc_time;        
        FSg_OS_dKO = sim_model_dKO.response_features(8).Overshoot;
        FSg_US_dKO = sim_model_dKO.response_features(8).Undershoot;
        FSg_npeak_dKO = length(sim_model_dKO.response_features(8).pks);
        FSg_Peak_dKO = sim_model_dKO.response_features(8).Peak;
        FSg_trise_dKO = sim_model_dKO.response_features(8).RiseTime;
        FSg_transient_dKO = sim_model_dKO.response_features(8).TransientTime;
        FSg_setmax_dKO = sim_model_dKO.response_features(8).SettlingMax;
        FSg_setmin_dKO = sim_model_dKO.response_features(8).SettlingMin;
        FSg_set_dKO = sim_model_dKO.response_features(8).SettlingTime;
        FSg_peaktime_dKO = sim_model_dKO.response_features(8).PeakTime;

        pMadcn_OScalc_dKO = sim_model_dKO.response_features(14).OS_calc;
        pMadcn_UScalc_dKO = sim_model_dKO.response_features(14).US_calc;
        pMadcn_OScalctime_dKO = sim_model_dKO.response_features(14).OS_calc_time;
        pMadcn_UScalctime_dKO = sim_model_dKO.response_features(14).US_calc_time;  
        pMadcn_OS_dKO = sim_model_dKO.response_features(14).Overshoot;
        pMadcn_US_dKO = sim_model_dKO.response_features(14).Undershoot;
        pMadcn_npeak_dKO = length(sim_model_dKO.response_features(14).pks);
        pMadcn_Peak_dKO = sim_model_dKO.response_features(14).Peak;
        pMadcn_trise_dKO = sim_model_dKO.response_features(14).RiseTime;
        pMadcn_transient_dKO = sim_model_dKO.response_features(14).TransientTime;
        pMadcn_setmax_dKO = sim_model_dKO.response_features(14).SettlingMax;
        pMadcn_setmin_dKO = sim_model_dKO.response_features(14).SettlingMin;
        pMadcn_set_dKO = sim_model_dKO.response_features(14).SettlingTime;

        Dadc_OScalc_dKO = sim_model_dKO.response_features(15).OS_calc;
        Dadc_UScalc_dKO = sim_model_dKO.response_features(15).US_calc;
        Dadc_OScalctime_dKO = sim_model_dKO.response_features(15).OS_calc_time;
        Dadc_UScalctime_dKO = sim_model_dKO.response_features(15).US_calc_time; 
        Dadc_OS_dKO = sim_model_dKO.response_features(15).Overshoot;
        Dadc_US_dKO = sim_model_dKO.response_features(15).Undershoot;
        Dadc_npeak_dKO = length(sim_model_dKO.response_features(15).pks);
        Dadc_Peak_dKO = sim_model_dKO.response_features(15).Peak;
        Dadc_trise_dKO = sim_model_dKO.response_features(15).RiseTime;
        Dadc_transient_dKO = sim_model_dKO.response_features(15).TransientTime;
        Dadc_setmax_dKO = sim_model_dKO.response_features(15).SettlingMax;
        Dadc_setmin_dKO = sim_model_dKO.response_features(15).SettlingMin;
        Dadc_set_dKO = sim_model_dKO.response_features(15).SettlingTime;
        Dadc_peaktime_dKO = sim_model_dKO.response_features(15).PeakTime; 

        FSc_OScalc_dKO = sim_model_dKO.response_features(16).OS_calc;
        FSc_UScalc_dKO = sim_model_dKO.response_features(16).US_calc;
        FSc_OScalctime_dKO = sim_model_dKO.response_features(16).OS_calc_time;
        FSc_UScalctime_dKO = sim_model_dKO.response_features(16).US_calc_time; 
        FSc_OS_dKO = sim_model_dKO.response_features(16).Overshoot;
        FSc_US_dKO = sim_model_dKO.response_features(16).Undershoot;
        FSc_npeak_dKO = length(sim_model_dKO.response_features(16).pks);
        FSc_Peak_dKO = sim_model_dKO.response_features(16).Peak;
        FSc_trise_dKO = sim_model_dKO.response_features(16).RiseTime;
        FSc_transient_dKO = sim_model_dKO.response_features(16).TransientTime;
        FSc_setmax_dKO = sim_model_dKO.response_features(16).SettlingMax;
        FSc_setmin_dKO = sim_model_dKO.response_features(16).SettlingMin;
        FSc_set_dKO = sim_model_dKO.response_features(16).SettlingTime;
        FSc_peaktime_dKO = sim_model_dKO.response_features(16).PeakTime;        


%         gsc_dKO = sim_model_dKO.gsc_profile;
%         precb_dKO = sim_model_dKO.precb_profile;
        if yesplot
            ax11.YLim(2) = 1.5*max([sim_model_wt.pMadgn sim_model_wt.pMadcn...
                sim_model_dKO.pMadgn sim_model_dKO.pMadcn...
                sim_model_wt.Dadg sim_model_wt.Dadc...
                sim_model_dKO.Dadg sim_model_dKO.Dadc...
                sim_model_wt.FSg sim_model_wt.FSc...
                sim_model_dKO.FSg sim_model_dKO.FSc...
                ]);
            ax12.YLim(2) = 1.5*max([sim_model_wt.pMadgn sim_model_wt.pMadcn...
                sim_model_dKO.pMadgn sim_model_dKO.pMadcn...
                sim_model_wt.Dadg sim_model_wt.Dadc...
                sim_model_dKO.Dadg sim_model_dKO.Dadc...
                sim_model_wt.FSg sim_model_wt.FSc...
                sim_model_dKO.FSg sim_model_dKO.FSc...
                ]);
            title(t1,['profile# ',num2str(loc)])
            saveas(gcf,['Plots',filesep,'profile_',num2str(loc),'.png'])
            % close(f1)
        end
        %
        %Simulate with a Dppg value of 1+1e-2 to calc sensitivity
        %coefficient
        %
        Dppg = 1; Dppg_delta = 1+1e-2; %these numbers are only used to calc phi
        %to calc sensitivity at a differernt value of Doog/Dppg_delta
        %change it in simulate_model.m

        sim_model_wt_delta = simulate_model(sets,Dppc_wt,wt,[],[],delta);

%         pMadg_wt_delta = sim_model_wt_delta.pMadg;
        pMadgn_wt_delta = sim_model_wt_delta.pMadgn;
%         Dadg_wt_delta = sim_model_wt_delta.Dadg;
%         FSg_wt_delta = sim_model_wt_delta.FSg;        

%         pMadc_wt_delta = sim_model_wt_delta.pMadc;
        pMadcn_wt_delta = sim_model_wt_delta.pMadcn;
%         Dadc_wt_delta = sim_model_wt_delta.Dadc;
%         FSc_wt_delta = sim_model_wt_delta.FSc;

%         gsc_wt_delta = sim_model_wt_delta.gsc_profile;
%         precb_wt_delta = sim_model_wt_delta.precb_profile;

        phi_pMadn_gsc_wt = (Dppg./pMadgn_wt).*((pMadgn_wt_delta-pMadgn_wt)./(Dppg_delta-Dppg));
        phi_pMadcn_precb_wt = (Dppg./pMadcn_wt).*((pMadcn_wt_delta-pMadcn_wt)./(Dppg_delta-Dppg));

        %calc phi for dad^KO

        sim_model_dKO_delta = simulate_model(sets,Dppc_wt,isdadko,[],[],delta);
        
%         pMadg_dKO_delta = sim_model_dKO_delta.pMadg;
        pMadgn_dKO_delta = sim_model_dKO_delta.pMadgn;

%         pMadc_dKO_delta = sim_model_dKO_delta.pMadc;
        pMadcn_dKO_delta = sim_model_dKO_delta.pMadcn;
%         Dadc_dKO_delta = sim_model_dKO_delta.Dadc;
%         FSc_dKO_delta = sim_model_dKO_delta.FSc;

%         gsc_dKO_delta = sim_model_dKO_delta.gsc_profile;
%         precb_dKO_delta = sim_model_dKO_delta.precb_profile;

        phi_pMadn_gsc_dKO = (Dppg./pMadgn_dKO).*((pMadgn_dKO_delta-pMadgn_dKO)./(Dppg_delta-Dppg));
        phi_pMadcn_precb_dKO = (Dppg./pMadcn_dKO).*((pMadcn_dKO_delta-pMadcn_dKO)./(Dppg_delta-Dppg));

        if plotturning && (isturning_dKO && isturning_wt)
            plot_turning(Xwt,XdKO);
            title(['#',num2str(loc)]);
            nexttile(6)
            hold on
            plot(Dppc_wt,pMadcn_wt,'Marker','o','Color','k','MarkerSize',12,'LineWidth',2)
            hold on
            plot(Dppc_wt,pMadcn_dKO,'Marker','o','Color','b','MarkerSize',12,'LineWidth',2)
            saveas(gcf,['Plots',filesep,'bifurcation_',num2str(loc),'.png'])
            % close
        end

        %
        % check if at steady state Dad is lower in CB than pMadgc
        % (Remember this is the only place we use pMadgc and not pMadgn)
        %
        pMadg_higher = (pMadg_wt > Dadg_wt) & (Dadg_wt > FSg_wt);
        Dad_low = Dadc_wt < 0.05;
        FS_turns_on = FSc_wt > 1.5*FSg_wt;
%         Dad_higher = Dadg_wt > 1.2*FSg_wt;
        FSc_max = (FSc_wt > Dadc_wt) & (FSc_wt > pMadc_wt);

        %
        % check if pMadgn in dad^KO is on the higher branch
        %
        if isturning_wt
        wt_lower_branch = pMadcn_wt < pMadn_LPL_wt;
        wt_upper_branch = pMadcn_wt > pMadn_LPU_wt;

        dadKO_lower_branch = pMadcn_dKO < pMadn_LPL_dKO;
        dadKO_upper_branch = pMadcn_dKO > pMadn_LPU_dKO;

        higher_branch = wt_lower_branch & dadKO_upper_branch;
        else 
        wt_lower_branch = nan;
        wt_upper_branch = nan;

        dadKO_lower_branch = nan;
        dadKO_upper_branch = nan;

        higher_branch = nan;
        end
%         branch = [wt_lower_branch;wt_upper_branch;dadKO_lower_branch;dadKO_upper_branch];
        %
        % check if the parameter set is perfect
        %

        perfectness_score = sum([(k_in_out_pMad > 1.5),isturning_wt,isturning_dKO,(Dppc_wt > 0.2),pMadg_higher,Dad_low,FS_turns_on,higher_branch,FSc_max]);

         datetime = {datestr(now,'yyyy-mm-dd-HH-MM-SS')};
         runtime = toc;
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
            pMadn_ratio_wt,FS_ratio_wt,pMadn_ratio_dKO,FS_ratio_dKO,datetime,runtime),...
            filename,"WriteVariableNames",false,"WriteRowNames",true,"WriteMode","append");

end

save('pscreen_rep_pss.mat','pscreen','-v7.3')



