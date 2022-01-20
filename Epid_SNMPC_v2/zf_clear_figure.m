%%
%  File: mf_clear_figure.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
% 
function [fig] = zf_clear_figure(fignr, figname)
arguments
    fignr = 8080;
    figname = "Figure";
end

fig = figure(fignr);
fig.Name = figname;
delete(fig.Children);

end