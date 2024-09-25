function dcdt = dpp_gsc_motif(t,c,params,Dppg)
%unpack state vars
Tkvg = c(1);
Tkvag = c(2);
Madg = c(3);
Madgn = c(4);
pMadg = c(5);
pMadgn = c(6);
Dadg = c(7);
FSg = c(8);

%unpack params
u = params.u;
mu = params.mu;
K = params.KHill;
tau = params.tau;
mu_Tkv = mu(1);
mu_pMad = mu(2);

%write hill here
pMadK2g = (pMadgn/K(1)).^4/(1+(pMadgn/K(1)).^4); %transcriptional actn and translation of Dad 
pMadK4g = ((pMadg/K(2)).^4)./(1+(pMadg/K(2)).^4); 

%compartment eqns
%GSC
dTkvg = ((u(2) + u(3)) + u(1)*Tkvag - Dppg.*Tkvg  - u(2)*Tkvg)/tau(1); 
dTkvag = (Dppg.*Tkvg - u(1)*Tkvag - u(9)*Dadg.* Tkvag - u(12)* FSg .* Tkvag - u(3)*Tkvag)/tau(1);
dMadg = (mu_pMad - mu_Tkv*Madg.* Tkvag + u(7)*Madgn - u(8)*Madg - mu_pMad*Madg)/tau(2);
dMadgn = (u(6)*pMadgn + u(8)*Madg - u(7)*Madgn)/tau(2);
dpMadg = (mu_Tkv*Madg*Tkvag + u(5)*pMadgn - u(4)*pMadg - u(10)*Dadg*pMadg - u(11)*pMadg*FSg - mu_pMad*pMadg)/tau(2); 
dpMadgn = (u(4)*pMadg - u(5)*pMadgn - u(6)*pMadgn)/tau(2);
dDadg = (pMadK2g - Dadg)/tau(3);
dFSg = (1 - u(13)*pMadK4g*FSg - FSg)/tau(4);

dcdt = [dTkvg;dTkvag;dMadg;dMadgn;dpMadg;dpMadgn;dDadg;dFSg];
end