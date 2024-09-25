function J = calc_jacobian(ss,params,Dppg)

F = dpp_gsc_motif([],ss,params,Dppg);

%pre-allocate for performance

Fh = NaN(size(F,1),size(F,1));

h = 0.0001;
for i = 1:size(ss,1)

    ss(i) = ss(i) + h;
    Fh(:,i) = dpp_gsc_motif([],ss,params,Dppg);
    ss(i) = ss(i) - h;

end
J = (Fh - F)./h;
end