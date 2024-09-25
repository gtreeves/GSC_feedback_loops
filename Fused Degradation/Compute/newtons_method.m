function [x, fval] = newtons_method(a,params,Dppg)

%define newton step params
maxEvals = 5000;
step_tolerance = 10^-10;

%define variables
syms t
% syms Dpp(t)
syms x [1 8]

% Dpp(t) = params.Dppg;
F(x) = dpp_gsc_motif(t,x,params,Dppg);

%calculate Jacobian
J(x) = jacobian(F,x);

%guess values
a = [a(1);a(2);a(3);a(4);a(5);a(6);a(7);a(8)];



%initialize newton counts
i = 1;

Fsum = norm(double((F(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8)))));
while Fsum > step_tolerance
%         F_old = double(F(a(1),a(2),a(3)));

    x_new = a - J(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8))\F(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8));
    F_new = double(F(x_new(1),x_new(2),x_new(3),x_new(4),x_new(5),x_new(6),x_new(7),x_new(8)));
    a = double(x_new);

    Fsum = norm(double(F_new));


    i = i + 1;
    if i==maxEvals
        break
    end

end

x = double(a);
fval = double(Fsum);
    

end