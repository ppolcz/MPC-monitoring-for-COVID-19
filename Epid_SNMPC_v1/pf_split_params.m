function [Np,pi,gamma,rhoI,rhoA,alpha,eta,h,mu,delta,vsp] = pf_split_params(params,args)
arguments
    params
    args.ParamType {mustBeMember(args.ParamType,["num","sym","nom","exp"])};
end
%%
%  File: mf_split_params.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%  Reviewed on 2021. November 17. (2021b)
%

switch args.ParamType
    case "nom"
        fn = 'nom';
    case {"num","exp"}
        fn = 'exp';
    case "sym"
        fn = 'val';
end

Np    = params.Np.(fn);
pi    = params.pi.(fn);
gamma = params.gamma.(fn);
rhoI  = params.rhoI.(fn);
rhoA  = params.rhoA.(fn);
alpha = params.alpha.(fn);
eta   = params.eta.(fn);
h     = params.h.(fn);
mu    = params.mu.(fn);
delta = params.delta.(fn);
vsp   = params.vsp.(fn);

switch nargout
    case 0
        Np = [ pi gamma rhoI rhoA alpha eta h mu delta vsp ];
    case 1
        Np = [ Np pi gamma rhoI rhoA alpha eta h mu delta vsp ];
    case 2
        pi = [ pi gamma rhoI rhoA alpha eta h mu delta vsp ];
end

end
