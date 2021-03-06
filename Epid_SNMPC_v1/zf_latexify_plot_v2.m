function [ret] = zf_latexify_plot_v2(r,args,xargs)
arguments
    r
    args
    xargs.Show_XLabel = true;
    xargs.Show_Months_Bars = true;
    xargs.Show_Months = true;
end
%%
%  File: mf_latexify_plot_v1.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%

ax = gca;
YLim = ax.YLim;

grid on
box on
hold on
ylabel(args.ylabel,'Interpreter','latex','FontSize',12);
if xargs.Show_XLabel
    xlabel(r.t_AxisLabel,'Interpreter','latex','FontSize',12);
end
xlim(r.d_Lims)

Logger.latexify_axis(gca,11);
% xticks(r.t_Ticks);
% xticklabels(r.t_Labels);

if isempty(ax.Legend)
    legend
end

if xargs.Show_Months_Bars
    Line0_XData =  r.d_Ticks' + [0;0;0];
    Line0_YData = [YLim NaN]' + zeros(size(r.d_Ticks))';
    Line0 = plot(Line0_XData(:),Line0_YData(:),':','Color',[1 1 1]*0.85);
    if xargs.Show_Months
        ax.Legend.String{end} = 'months';
    end
end

Removable_Leg_Entries = strcmp(ax.Legend.String,'') | startsWith(ax.Legend.String,'data');
if all(Removable_Leg_Entries)
    delete(ax.Legend);
else
    ax.Legend.String(Removable_Leg_Entries) = [];
end

if ~isempty(ax.Legend)
    ax.Legend.Interpreter = 'latex';
    ax.Legend.FontSize = r.Leg_FontSize;
    % ax.Legend.Location = r.Leg_Location;
end

ax.YLim = [ max(0,YLim(1)) YLim(2) ];
ax.YTickMode = "auto";

drawnow

end