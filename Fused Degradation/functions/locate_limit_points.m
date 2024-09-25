function limit_points = locate_limit_points(X)



%there are two limit points: left and right, right is lower than left, and
%right is the first time lambda_dot turns

dX_dpp = diff(X(end,:),1);
i = 1;
try
    while dX_dpp(:,i) > 0
        i=i+1;
    end
    lowerloc = i-1;
    clear i
    i = lowerloc+1;
    while dX_dpp(:,i) < 0
        i=i+1;
    end
    upperloc = i + 1;
    lp_lower = X(end,lowerloc);
    lp_upper = X(end,upperloc);

    lp_avg = (lp_lower + lp_upper)/2;

    Dpp_LPL = X(end,lowerloc);
    Dpp_LPU = X(end,upperloc);
    pMad_LPL = X(6,lowerloc);
    pMad_LPU = X(6,upperloc);
    %

    limit_points.Dpp_LPL = Dpp_LPL;
    limit_points.Dpp_LPU = Dpp_LPU;
    limit_points.pMad_LPL = pMad_LPL;
    limit_points.pMad_LPU = pMad_LPU;
end
end


