function r = pf_resolve_params(Par,params,r)
arguments
    Par
    params
    r = struct;
end
%%
%  File: mf_resolve_params.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%

nom = 'nom';          % kozepertek.
pm_perc = 'pm_perc';  % [1-pm_perc,1+pm_perc]*nom 0.95%-os szig. int.
pm_val = 'pm_val';    % [-pm_val,+pm_val]+nom 0.95%-os szig. int.
lim = 'lim';          % 0.95%-os szignifikancia intervallum (2 szigma)


fns = fieldnames(Par);

for i = 1:numel(fns)
    fn = fns{i};

    Par.(fn).unc = false;
    Par.(fn).var = 0;
    if isfield(Par.(fn),nom)

        Par.(fn).exp = Par.(fn).nom;
        if isfield(Par.(fn),lim)
            Par.(fn).exp = mean(Par.(fn).lim);
            Par.(fn).pm_val = Par.(fn).exp - min(Par.(fn).lim);
        end        

        if isfield(Par.(fn),pm_perc)
            Par.(fn).pm_val = Par.(fn).exp * Par.(fn).pm_perc / 100;
        end
        
        if isfield(Par.(fn),pm_val)
            Par.(fn).unc = true;
            Par.(fn).var = ( Par.(fn).pm_val / 2 )^2;
            Par.(fn).pm_perc = Par.(fn).pm_val / Par.(fn).exp * 100;
            Par.(fn).lim = [-1 1]*Par.(fn).pm_val + Par.(fn).exp;
        end

    else
        Par.(fn).nom = NaN;
        Par.(fn).exp = NaN;
        Par.(fn).lim = [NaN NaN];
        Par.(fn).pm_perc = NaN;
        Par.(fn).pm_val = NaN;
    end

    Par.(fn).val = Par.(fn).exp;
end

np = numel(params);
Mp = zeros(np,1);
Sp = zeros(np,1);
Up = ~Sp;
for i = 1:np
    s = params(i);
    fn = s.name;
    p = Par.(fn);
    Mp(i) = p.exp;
    Sp(i) = p.var;
    Up(i) = p.unc;

    if Up(i)
        Par.(fn).val = s;
    end
end


% Updated parameter information
r.Par = Par;

% Number of uncertain parameters
r.np = numel(params);

% Mean value of the parameters: mu_theta (theta == params)
r.Mp = Mp(Up);

% Variance matrix of the parameters: Sigma_theta
r.Sp = diag(Sp(Up));

% Select out deterministic parameters. (Use find: CasADi does not provide
% logical indexing.)
r.params = params(find(Up)); %#ok<FNDSB> 

end