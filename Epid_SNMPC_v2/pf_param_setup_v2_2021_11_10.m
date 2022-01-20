function r = pf_param_setup_v2_2021_11_10(r)
arguments
    r = struct;
end
%%
%  File: pf_param_setup_v2_2021_11_10.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%

import casadi.*

% Population in Hungary
Par.Np.nom = 9800000;

% Inverse of pre-symptomatic infectious period [1/days]
Par.pi.nom = 1/3;
% Par.pi.lim = [1/2,1/4];
Par.pi.pm_perc = 30;

% Probability of developing symptoms
Par.gamma.alt_values = 0.75;
Par.gamma.nom = 0.6;
Par.gamma.pm_perc = 10;

% Inverse of infectious period [1/days]
Par.rhoI.nom = 1/4;
% Par.rhoI.lim = [1/3,1/5];
Par.rhoI.pm_perc = 25;

% Inverse of infectious period [1/days]
Par.rhoA.nom = 1/4;
% Par.rhoA.lim = [1/3,1/5];
Par.rhoA.pm_perc = 25;


% Inverse of latent period [1/days]
Par.alpha.nom = 1/2.5;
% Par.alpha.lim = [1/2,1/3];
Par.alpha.pm_perc = 20;

% Hospitalization probability of symptomatic cases
% (Feltételezhető, hogy ezt elég pontosan tudjunk)
Par.eta.nom = 0.076;
Par.eta.pm_perc = 10;

% Inverse of average hospitalization length [1/days]
% (Feltételezhető, hogy ezt elég pontosan tudjunk)
Par.h.nom = 1/10;
% Par.h.lim = [1/9,1/11];
Par.h.pm_perc = 10;

% Probability of fatal outcome in hospital
Par.mu.alt_values = [0.145 , 0.185];
Par.mu.nom = 0.205; % 2021.09.28. (szeptember 28, kedd), 11:16
Par.mu.pm_perc = 10;

% Relative transmissibility of asymptomatic
Par.delta.nom = 0.75;
Par.delta.pm_perc = 10;

% Vaccination success percentage
Par.vsp.nom = 0.75;
Par.vsp.pm_perc = 10;

Par.beta.var_Alpha = 1/3;
Par.beta.var_Britt = 0.5;
Par.beta.var_Delta = 0.7;

params = [
    SX.sym('pi')
    SX.sym('gamma')
    SX.sym('rhoI')
    SX.sym('rhoA')
    SX.sym('alpha')
    SX.sym('eta')
    SX.sym('h')
    SX.sym('mu')
    SX.sym('delta')
    SX.sym('vsp')
    ];

r = pf_resolve_params(Par,params,r);

end