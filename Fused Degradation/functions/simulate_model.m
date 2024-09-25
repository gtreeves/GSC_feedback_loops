function sim_model = simulate_model(sets,lp_avg,isdadko,yesplot,newwindow,delta)

%
%if dadko not provided
%
switch nargin
    case 2
        isdadko = 0;
        yesplot = 0;
        newwindow = 0;
        delta = 0;
    case 3
        yesplot = 0;
        newwindow = 0;
        delta = 0;
    case 4
        newwindow = 0;
        delta = 0;
    case 5 
        delta = 0;
end
%
% if isdadko, make parameter#16 really high
%
%

if isdadko
    sets(:,16) = 10000;
end

if delta
    Dppg = 1+1e-2;
else
    Dppg = 1;
end
%
% unpack sets
%
params.u = sets(:,1:13);
params.mu = sets(:,14:15);
params.KHill = sets(:,16:17);
% params.lambda = sets(:,19);
% params.alpha = sets(:,20);
params.tau = sets(:,18:21);
%
% Stuff required to run dpp_gsc_motif
%
fhandle = @dpp_gsc_motif;
c0 = zeros(1,8);
tspan = [0 43200]; %12*60*60
tauVF = 1000;
% Dppg = 1;
Dppc = lp_avg;
km = @(t)0.03*(exp(-t/tauVF));
km0 = @(t)0.03;

% Simulate GSC to steady state
[t,c] = ode15s(fhandle,tspan,c0,[],params,Dppg);
steadyC = newtons_method(c(end,:),params,Dppg);

% Simulate for another hour with mass transfer (the boundary does not
% shrink during this interval, the cell grows)
[tc0,cc0] = ode15s(@dpp_gsc_compartments,2*tspan,[steadyC,steadyC],[],params,Dppg,Dppc,km0);

% Simulate with mass transfer and boundary shrinks acc. to tauVF
[tc,cc] = ode15s(@dpp_gsc_compartments,tspan,cc0(end,:),[],params,Dppg,Dppc,km);

pMadg = cc(end,5);
pMadgn = cc(end,6);
Dadg = cc(end,7);
FSg = cc(end,8);

pMadc = cc(end,13);
pMadcn = cc(end,14);
Dadc = cc(end,15);
FSc = cc(end,16);

if yesplot
    if newwindow
        figure
    end
    plot(tc/3600,cc(:,6),'b-',tc/3600,cc(:,14),'b--',linewidth=2)
    hold on
    plot(tc/3600,cc(:,8),'g-',tc/3600,cc(:,16),'g--',linewidth=2)
    hold on
    plot(tc/3600,cc(:,7),'r-',tc/3600,cc(:,15),'r--',linewidth=2)
    hold on
    xline(1.5)
    hold off
    grid on
    grid minor
    
    lgd = legend('pMad','','FS','','Dad','','Location','northoutside','NumColumns',3);
    title(lgd,['— GSC',' -- preCB'])
    xlabel('Time (hr)')
    ylabel('Concentration')
    set(gca,'FontSize',14)
    ylim([0 max(cc(:,[6,7,8,14,15,16]),[],'all')])
    xlim([0 tc(end)/3600])
end

response_features = characterize_signal(cc,tc);

sim_model.pMadg = pMadg;
sim_model.pMadgn = pMadgn;
sim_model.Dadg = Dadg;
sim_model.FSg = FSg;

sim_model.pMadc = pMadc;
sim_model.pMadcn = pMadcn;
sim_model.Dadc = Dadc;
sim_model.FSc = FSc;

sim_model.gsc_profile = [c,t];
sim_model.mitosis = [cc0,tc0];
sim_model.precb_profile = [cc,tc];
sim_model.response_features = response_features;
end

