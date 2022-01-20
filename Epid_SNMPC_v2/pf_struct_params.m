function s = pf_struct_params(params,args)
arguments
    params
    args.ParamType {mustBeMember(args.ParamType,["num","sym","nom","exp"])};
end
%%
%  File: pf_struct_params.m
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

s.Np    = params.Np.(fn);
s.pi    = params.pi.(fn);
s.gamma = params.gamma.(fn);
s.rhoI  = params.rhoI.(fn);
s.rhoA  = params.rhoA.(fn);
s.alpha = params.alpha.(fn);
s.eta   = params.eta.(fn);
s.h     = params.h.(fn);
s.mu    = params.mu.(fn);
s.delta = params.delta.(fn);
s.vsp   = params.vsp.(fn);

end
