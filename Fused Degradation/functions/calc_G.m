function G = calc_G(ss,params,Dppg)

F = dpp_gsc_motif([],ss,params,Dppg);
h = 0.0001;
Gh = dpp_gsc_motif([],ss,params,Dppg+h);
G = (Gh-F)./h;

end