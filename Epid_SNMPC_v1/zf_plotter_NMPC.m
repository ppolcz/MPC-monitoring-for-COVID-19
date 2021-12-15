function zf_plotter_NMPC(r,args)
arguments
    r
    args.Layout = {2 3};
    args.FigSize = [1900 980];
    args.FigNr = 1;
    args.FigName = 'Figure';
    args.Version = '0';
    args.Plotter = @zf_latexify_plot_v2;
end
%%
%  File: zf_plotter_NMPC.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study3_Gabor_kodjai_val
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2021. February 17. (2020b)
%  Minor review on 2021. September 15. (2021a)
%  Reviewed on 2021. October 11. (2021b)
%  Reviewed on 2021. October 12. (2021b)
%

%%

FilePref = sprintf(['results/NMPC%s/%s_Pred-%s'],...
    args.Version,datestr(date,29),datestr(r.Date_End_REC,29));
DIR = fileparts(FilePref);
if ~exist(DIR,"dir")
    mkdir(DIR)
end

r = zf_plotter_opts(r);

N = r.N_rec;

%% REC: Segedvaltozok

L_val = r.x_val(2,:)';
P_val = r.x_val(3,:)';
I_val = r.x_val(4,:)';
A_val = r.x_val(5,:)';
H_val = r.x_val(6,:)';
R_val = r.x_val(7,:)';
D_val = r.x_val(8,:)';
U_val = r.x_val(9,:)';
R_val_val = r.x_val(10,:)';

%% Allitsuk be a Figure ablakot floating mode-ba.

fig = zf_clear_figure(args.FigNr,args.FigName);
fig.Position([3,4]) = args.FigSize;
Tl = tiledlayout(args.Layout{:});
Tl.TileSpacing = "compact";
Tl.Padding = "compact";
Tl_N = prod(Tl.GridSize);

%%

Perc_Max = 0.25;

Idx = 4;
if Idx <= Tl_N
    Ax(Idx) = nexttile(Idx);
        
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',14);
    ylim([0,Perc_Max])
    
    yyaxis left
    plot(r.d_val,r.x_val(6,:),'.-','Color',[0.4660 0.6740 0.1880]), hold on
    plot(r.d_val(1:N+1),r.H_ref_All(1:N+1)','-r','LineWidth',1.5)
    legend('($\mathbf{H}$) Closed-loop trajectory for fixed parameters','($\mathbf{H}^{\mathrm{Off}}$) Official hospitalization data','Location','northwest')
    
    title('($\mathbf{H}$) Hosp. compared to official data','Interpreter','latex','FontSize',14)
    
    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    args.ylabel = '';
    args.fname = '_Skip_H';
    args.Plotter(r,args,"Show_XLabel",Idx > (args.Layout{1}-1)*args.Layout{2})
end

%%

Perc_Max = 100;

Idx = 3;
if Idx <= Tl_N
    Ax(Idx) = nexttile(Idx);
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',14);
    ylim([0,Perc_Max])
    
    yyaxis left
    plot(r.d_val,R_val,'LineWidth',2);
    hold on
    plot(r.d_val,R_val_val,'LineWidth',2);
    plot(r.d_val,R_val + U_val,'LineWidth',2);
    plot(0,0,'w')
    plot(0,0,'w')
    
    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title('Number of immune people','Interpreter','latex','FontSize',14)
    
    legend(...
        'Immune by recovery ONLY',...
        'Total recovereds',...
        sprintf('Immune by recovery and/or vaccination'),...
    ... sprintf('(prediction from %s: %d/day)',datestr(t_date_REC(Idx_Last_Available_Data+Oltasok_kozotti_ido),'mmm. dd, yyyy'),V_(end)),...
        'Location','northwest');
    
    args.ylabel = 'no. of immune people';
    args.fname = '_Fig3_recovered';
    args.Plotter(r,args,"Show_XLabel",Idx > (args.Layout{1}-1)*args.Layout{2})
end

%%

Idx = 2;
if Idx <= Tl_N
    Ax(Idx) = nexttile(Idx);
    
    Pb_beta = plot(r.d_val(1:N),r.u_val(1:N)); hold on
    Pb_Rc = plot(r.d_val(1:N),r.Rc_t(1:N));
    
    ylim([0 3])
    
    title('($\beta$) Transmission rate (infectiousness)','Interpreter','latex','FontSize',14)
    
    legend([Pb_beta,Pb_Rc], {'($\beta$) Transmission rate','($R_t$) Time-dependent reproduction number'},...
        'Interpreter','latex','Location','northwest')
    
    args.ylabel = '($\beta$) Transmission rate';
    args.fname = '_Fig2_beta';
    args.Plotter(r,args,"Show_XLabel",Idx > (args.Layout{1}-1)*args.Layout{2})
end

%%

Perc_Max = 1.5;

Idx = 1;
if Idx <= Tl_N
    Ax(Idx) = nexttile(Idx);
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',14);
    ylim([0,Perc_Max])

    yyaxis left
    hold on
    plot(r.d_val,L_val,'-','Color',[0 0.4470 0.7410])
    plot(r.d_val,P_val,'-','Color',[0.8500 0.3250 0.0980])
    plot(r.d_val,I_val,'-','Color',[0.9290 0.6940 0.1250])
    plot(r.d_val,A_val,'-','Color',[0.4940 0.1840 0.5560])
    plot(r.d_val,H_val,'-','Color',[0.4660 0.6740 0.1880])
    plot(r.d_val,D_val,'-','Color',[0.3010 0.7450 0.9330])
    plot(r.d_val(1:N+1),r.H_ref_All(1:N+1)','-r','LineWidth',1.5)
    
    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title('($\mathbf{L}$--$\mathbf{H},\mathbf{D}$) Est. trajectories for fixed parameters',...
        'Interpreter','latex','FontSize',10)
    Leg = legend('$\textbf{L}$','$\textbf{P}$','$\textbf{I}$','$\textbf{A}$','$\textbf{H}$','$\textbf{D}$',...
        '$\mathbf{H}^{\mathrm{Off}}$ (Official hospitalization data)',...
        'Interpreter','latex','Location','northwest');
    
    Leg.NumColumns = 4;
    
    args.ylabel = '';
    args.fname = '_Fig1_LPIAHD';
    args.Plotter(r,args,"Show_XLabel",Idx > (args.Layout{1}-1)*args.Layout{2})
end

%%

Perc_Max = 5;

Idx = 5;
if Idx <= Tl_N
    Ax(Idx) = nexttile(Idx);

    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',14);
    ylim([0,Perc_Max])

    yyaxis left
    plot(r.d_val,r.Daily_all_Inf');

    legend('All infected',...
        'Interpreter','latex','Location','northwest')

    ylim([0,Perc_Max]*r.Par.Np.val/100);

    title('All infected ($\textbf{L} + \textbf{P} + \textbf{I} + \textbf{A} + \textbf{H}$)',...
        'Interpreter','latex','FontSize',11)

    args.ylabel = '';
    args.fname = '_Fig5_val_inf';
    args.Plotter(r,args,"Show_XLabel",Idx > (args.Layout{1}-1)*args.Layout{2})
end

%%

Perc_Max = 0.5;

Idx = 6;
if Idx <= Tl_N
    Ax(Idx) = nexttile(Idx);

    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',14);
    ylim([0,Perc_Max])
    
    yyaxis left
    plot(r.d_val,r.Daily_new_Inf);

    ylim([0,Perc_Max]*r.Par.Np.val/100);

    legend('Daily new infected',...
        'Interpreter','latex','Location','northwest')

    title(['Daily new infected ' ...
        '($\beta_k(\textbf{P}_k + \textbf{I}_k + \delta \textbf{A}_k) \textbf{S}_k / \textbf{N} $)'],...
        'Interpreter','latex','FontSize',11)

    args.ylabel = '';
    args.fname = '_Fig6_Daily_new';
    args.Plotter(r,args,"Show_XLabel",Idx > (args.Layout{1}-1)*args.Layout{2})
end

%%

OBJ_str = num2str(r.NMPC_sol.F.H_error_SUM);

% keyboard
pause(0.5);

exportgraphics(fig,[ FilePref '_Fig__' OBJ_str '.png' ])

% for i = 1:numel(Ax)
%     exportgraphics(Ax(i),[ FilePref '_Ax' num2str(i) '__' OBJ_str '.png' ])
% end

saveas(fig,[ FilePref '_Fig__' OBJ_str '.fig' ])

end
