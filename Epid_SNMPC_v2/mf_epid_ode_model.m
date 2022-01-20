function [x,beta,V,f,C,h,hvar] = mf_epid_ode_model(r,args)
arguments
    r
    args.ParamType {mustBeMember(args.ParamType,["num","sym","nom","exp"])};
end
%%
%  File: mf_epid_ode_model.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%

import casadi.*

[Np,pi,gamma,rhoI,rhoA,alpha,eta,h,mu,delta,vsp] = pf_split_params(r.Par,"ParamType",args.ParamType);

% State variables:
S = SX.sym('S'); % 1
L = SX.sym('L'); % 2
P = SX.sym('P'); % 3
I = SX.sym('I'); % 4
A = SX.sym('A'); % 5
H = SX.sym('H'); % 6
R = SX.sym('R'); % 7
D = SX.sym('D'); % 8
U = SX.sym('U'); % 9
% -----
R_All = SX.sym('R_All'); % 10
% -----
x = [S;L;P;I;A;H;R;D;U;R_All];
C = [0,0,0,0,0,1,0,0,0,0];

% Time-dependent parameter: Nr. of vaccinated per day
V = SX.sym('V');

% Control input:
beta = SX.sym('beta');

% 2021.09.23. (szeptember 23, csütörtök), 00:04         |
V_mlt = vsp * V / (Np - H - I - D - U - L - P - A); %   |
                                        %               V
dS = -beta * (P + I + delta*A) * S / Np                 - S * V_mlt;
dL =  beta * (P + I + delta*A) * S / Np - alpha * L;
dP = alpha * L - pi * P;
dI = gamma * pi * P - rhoI * I;
dA = (1-gamma) * pi * P - rhoA * A;
dH = rhoI * eta * I - h * H;
dR = rhoI * (1 - eta) * I + rhoA * A + (1 - mu) * h * H  - R * V_mlt;
dD = mu * h * H;
dU = vsp * V;
% -----
dR_All = rhoI * (1 - eta) * I + rhoA * A + (1 - mu) * h * H;

dxdt = [dS dL dP dI dA dH dR dD dU dR_All].';

Ts = 1;
f = x + Ts*dxdt;

h = {};
h.val = [
    L + P + I + A + H
    beta * (P + I + delta*A) * S / Np
    beta * (1/pi + gamma/rhoI + delta*(1-gamma)/rhoA) * S / Np
    ];

Is_Stochastic = strcmp(args.ParamType,"sym");

if Is_Stochastic
    h.Fn = Function('h',{x,beta,r.params},{h.val},{'x','u','theta'},{'[all_inf,daily_new_inf,Rc]'});
else
    h.Fn = Function('h',{x,beta},{h.val},{'x','u'},{'[all_inf,daily_new_inf,Rc]'});
end

h.desc_Input_1 = 'State vector (x)';
h.desc_Input_2 = 'Input vector (u = beta)';
if Is_Stochastic
    h.desc_Input_3 = 'Parameter vector (theta)';
end

h.desc_Output_dim1 = 'All infected';
h.desc_Output_dim2 = 'Daily new infected';
h.desc_Output_dim3 = 'Time-dependent reproduction number (Rc)';

hvar = {};
hvar.val = zeros(size(h.val));
hvar.Fn = [];
if Is_Stochastic
    vars = [x;beta;r.params];
    Jh = jacobian(h.val,vars);
    
    [Sigma,Sigma_half] = Pcz_CasADi_Helper.create_sym('Sigma',numel(vars),1,'str','sym');

    hvar.val = diag(Jh * Sigma * Jh');
    hvar.Fn = Function('hvar',{x,beta,r.params,Sigma_half},{hvar.val},{'x','u','theta','Sigma'},{'[Var_all_inf,Var_daily_new_inf,Var_Rc]'});
end
end