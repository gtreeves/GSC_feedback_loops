function data = run_palc(sets,isdadko)

switch nargin
    case 1
        isdadko = 0;
end
%PALC parameters
step_size = 0.01;
% maxEvalsnewton = 5000;
maxEvalspalc = 100000;
step_tolerance_newton = 10^-10;


%initial conditions for the model
c0 = zeros(1,8);
tspan = [0 5000]; %just used to give newton a good guess to find steady state, so not running for 12*60*60

Dppg = 0.01; %can be a time-dependant function
Dppmin = 0; %choose how long you want to run your code for
Dppmax = 1; %choose how long you want to run your code for

if isdadko
    sets(:,15) = 10000;
end

params.u = sets(:,1:12);
params.mu = sets(:,13:14);
params.KHill = sets(:,15:16);
% params.lambda = sets(:,19);
% params.alpha = sets(:,20);
params.tau = sets(:,17:20);

[~,c] = ode15s(@dpp_gsc_motif,tspan,c0,[],params,Dppg);

ss = newtons_method(c(end,:),params,Dppg);

%This is regular continuation
Fcalc = dpp_gsc_motif([],ss,params,Dppg);
Jcalc = calc_jacobian(ss,params,Dppg);
Gcalc = calc_G(ss,params,Dppg);

% for first predictor
dxdl = -Jcalc\Gcalc; % GTR: this has a negative sign
lambda_dot = 1/(sqrt(1+(norm(dxdl))^2)); %always positive
x_dot = dxdl*lambda_dot;

%guess for newton's
x_guess = ss + x_dot*step_size;
lambda_guess = Dppg + lambda_dot*step_size; %here
xp = [x_guess;lambda_guess];

xsoln = ss; % new
lambdasoln = Dppg; % new

%PALC starts here
%now start a palc counter and save these (x,lambda) values

palc_counter = 1;
L_DOT(palc_counter,:) = lambda_dot;
X(:,palc_counter) = [xsoln;lambdasoln];
palc_counter = palc_counter + 1; %update counter
%create an elementary vector
en = [zeros(size(x_guess,1),1);ones(1,1)];


while lambdasoln > Dppmin && lambdasoln < Dppmax
    
	Fsum = 1; % dummy value to enter the newton loop
    ii = 1;
    while Fsum > step_tolerance_newton

        %
		% Calc F
		%
        Fcalc = dpp_gsc_motif([],xp(1:end-1),params,xp(end));            
        F_hat = [Fcalc;(x_dot'*(xp(1:end-1)-xsoln)+lambda_dot*(xp(end)-lambdasoln) - step_size)]; %here
		
		%
		% Calc Jacobian
		%
        Jcalc = calc_jacobian(xp(1:end-1),params,xp(end));
        Gcalc = calc_G(xp(1:end-1),params,xp(end));
        J_hat = [Jcalc,Gcalc;x_dot',lambda_dot];
		
        x_new = xp - J_hat\F_hat;
        xp = x_new;
        
        Fsum = norm(F_hat);

%         ii = ii + 1;
%         if ii==maxEvalsnewton
%             disp('Maximum Evaluations Exceeded')
%             break
%         end

    end
    clear ii
    xsoln = xp(1:end-1);
    lambdasoln = xp(end);
    
    bigX = J_hat\en;
    bigX = bigX./norm(bigX);
    x_dot = bigX(1:end-1);
    lambda_dot = bigX(end);

    xp = xp + bigX*step_size;

    L_DOT(palc_counter,:) = lambda_dot;
    X(:,palc_counter) = xp;
    palc_counter = palc_counter + 1;
    
    if any(L_DOT<0)
        isturning = 1;
    else
        isturning = 0;
        
    end
    
     if palc_counter==maxEvalspalc
         disp('PALC exceeded maxEvals')
         break
     end

    data.isturning = isturning;

%     if isturning 

data.X = X;
data.L_DOT = L_DOT;

%     else
%         
%     data.X = [];
%     data.L_DOT = [];
% 
%     end

end