function [ret] = zf_latexify_plot_v1(r,args)
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
ylabel(args.ylabel,'Interpreter','latex','FontSize',14);
xlabel(r.t_AxisLabel,'Interpreter','latex','FontSize',14);
xlim(r.d_Lims)

Logger.latexify_axis(gca,14);
% xticks(r.t_Ticks);
% xticklabels(r.t_Labels);

Line0_XData =  r.d_Ticks' + [0;0;0];
Line0_YData = [YLim NaN]' + zeros(size(r.d_Ticks))';
Line0 = plot(Line0_XData(:),Line0_YData(:),':','Color',[1 1 1]*0.85);

Line1 = line([0 0]+r.d_Start_MPC,ax.YLim,'Color','blue','LineStyle','--');
Line2 = line([0 0]+r.d_Last_Available_Data,ax.YLim,'Color','black','LineStyle','--','LineWidth',1);
Line3 = line([0 0]+r.d_Current_Date,ax.YLim,'Color','black','LineStyle','--','LineWidth',2);

Line1.YData = YLim;
Line2.YData = YLim;
Line3.YData = YLim;

if isempty(ax.Legend)
    legend
    ax.Legend.String{1} = args.ylabel;
end
ax.Legend.String{end-3} = 'months';
ax.Legend.String{end-2} = r.Str_Date_Start_MPC;
ax.Legend.String{end-1} = r.Str_Last_Available_Data;
ax.Legend.String{end-0} = r.Str_Current_Date;

if ~isempty(ax.Legend)
    ax.Legend.Interpreter = 'latex';
    ax.Legend.FontSize = r.Leg_FontSize;
end

ax.Legend.Location = r.Leg_Location;

ax.YLim = [ max(0,YLim(1)) YLim(2) ];
ax.YTickMode = "auto";

drawnow

end