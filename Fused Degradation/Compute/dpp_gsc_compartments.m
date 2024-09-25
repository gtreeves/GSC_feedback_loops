function dcdt = dpp_gsc_compartments(t,c,params,Dppg,Dppc,km)
%unpack state vars
Tkvg = c(1);
Tkvag = c(2);
Madg = c(3);
Madgn = c(4);
pMadg = c(5);
pMadgn = c(6);
Dadg = c(7);
FSg = c(8);

Tkvc = c(9);
Tkvac = c(10);
Madc = c(11);
Madcn = c(12);
pMadc = c(13);
pMadcn = c(14);
Dadc = c(15);
FSc = c(16);

%unpack params
u = params.u;
mu = params.mu;
K = params.KHill;
tau = params.tau;

mu_Tkv = mu(1);
mu_pMad = mu(2);

%write Hills here
%gsc
pMadK2g = (pMadgn/K(1)).^4/(1+(pMadgn/K(1)).^4); %transcriptional actn and translation of Dad 
pMadK4g = ((pMadg/K(2)).^4)./(1+(pMadg/K(2)).^4); %transcriptional actn and repressed translation of Fu-S
%cb
pMadK2c = (pMadcn/K(1)).^4/(1+(pMadcn/K(1)).^4); %transcriptional actn and translation of Dad 
pMadK4c = ((pMadc/K(2)).^4)./(1+(pMadc/K(2)).^4); %transcriptional actn and repressed translation of Fu-S

%compartment eqns
%GSC
dTkvg = ((u(2) + u(3)) + u(1)*Tkvag - Dppg.*Tkvg  - u(2)*Tkvg)/tau(1); 
dTkvag = (Dppg.*Tkvg - u(1)*Tkvag - u(9)*Dadg.* Tkvag - u(12)* FSg .* Tkvag - u(3)*Tkvag)/tau(1);
dMadg = (mu_pMad - mu_Tkv*Madg.* Tkvag + u(7)*Madgn - u(8)*Madg - mu_pMad*Madg)/tau(2) - km(t)*(Madg-Madc);
dMadgn = (u(6)*pMadgn + u(8)*Madg - u(7)*Madgn)/tau(2);
dpMadg = (mu_Tkv*Madg*Tkvag + u(5)*pMadgn - u(4)*pMadg - u(10)*Dadg*pMadg - u(11)*pMadg*FSg - mu_pMad*pMadg)/tau(2) - km(t)*(pMadg-pMadc); 
dpMadgn = (u(4)*pMadg - u(5)*pMadgn - u(6)*pMadgn)/tau(2);
dDadg = (pMadK2g - Dadg)/tau(3) - km(t)*(Dadg-Dadc);
dFSg = (1 - u(13)*pMadK4g*FSg - FSg)/tau(4) - km(t)*(FSg-FSc);

%CB
dTkvc = ((u(2) + u(3)) + u(1)*Tkvac - Dppc.*Tkvc  - u(2)*Tkvc)/tau(1); 
dTkvac = (Dppc.*Tkvc - u(1)*Tkvac - u(9)*Dadc.* Tkvac - u(12)* FSc .* Tkvac - u(3)*Tkvac)/tau(1);
dMadc = (mu_pMad - mu_Tkv*Madc.* Tkvac + u(7)*Madcn - u(8)*Madc - mu_pMad*Madc)/tau(2) + km(t)*(Madg-Madc);
dMadcn = (u(6)*pMadcn + u(8)*Madc - u(7)*Madcn)/tau(2);
dpMadc = (mu_Tkv*Madc*Tkvac + u(5)*pMadcn - u(4)*pMadc - u(10)*Dadc*pMadc - u(11)*pMadc*FSc - mu_pMad*pMadc)/tau(2) + km(t)*(pMadg-pMadc); 
dpMadcn = (u(4)*pMadc - u(5)*pMadcn - u(6)*pMadcn)/tau(2);
dDadc = (pMadK2c - Dadc)/tau(3) + km(t)*(Dadg-Dadc);
dFSc = (1 - u(13)*pMadK4c*FSc - FSc)/tau(4) + km(t)*(FSg-FSc);

dcdt = [dTkvg;dTkvag;dMadg;dMadgn;dpMadg;dpMadgn;dDadg;dFSg;dTkvc;dTkvac;dMadc;dMadcn;dpMadc;dpMadcn;dDadc;dFSc];
end