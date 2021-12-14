function zf_plotter_MDPI(r,args)
arguments
    r
    args.Layout = {3 2};
    args.FigSize = [1053 1256];
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

FilePref = sprintf(['results/MDPI%s/%s_Pred-%s'],...
    args.Version,datestr(date,29),datestr(r.Date_End_REC,29));
DIR = fileparts(FilePref);
if ~exist(DIR,"dir")
    mkdir(DIR)
end

r = zf_plotter_opts(r);

%% REC: Segedvaltozok

d_recpred = r.d_recpred';
N_rec = r.N_rec;
H_ref = r.H_ref;

Exp_u = r.Muu_val;
Std_u = sqrt(r.Suu_val);

Mxx_val = r.Mxx_val;
Exp_S = Mxx_val(1,:)';
Exp_L = Mxx_val(2,:)';
Exp_P = Mxx_val(3,:)';
Exp_I = Mxx_val(4,:)';
Exp_A = Mxx_val(5,:)';
Exp_H = Mxx_val(6,:)';
Exp_R = Mxx_val(7,:)';
Exp_D = Mxx_val(8,:)';
Exp_U = Mxx_val(9,:)';
Exp_R_All = Mxx_val(10,:)';

Sxx_diag = r.Sxx_diag;
Std_S = sqrt(Sxx_diag(1,:)');
Std_L = sqrt(Sxx_diag(2,:)');
Std_P = sqrt(Sxx_diag(3,:)');
Std_I = sqrt(Sxx_diag(4,:)');
Std_A = sqrt(Sxx_diag(5,:)');
Std_H = sqrt(Sxx_diag(6,:)');
Std_R = sqrt(Sxx_diag(7,:)');
Std_D = sqrt(Sxx_diag(8,:)');
Std_U = sqrt(Sxx_diag(9,:)');
Std_R_All = sqrt(Sxx_diag(10,:)');

%% Allitsuk be a Figure ablakot floating mode-ba.

fig = zf_clear_figure(args.FigNr,args.FigName);
fig.Position([3,4]) = args.FigSize;
Tl = tiledlayout(args.Layout{:});
Tl.TileSpacing = "compact";
Tl.Padding = "compact";

%% Titles

Cnt(0);
Titles = {
    [ '\textbf{Plot ' num2str(Cnt) '.} Expected values' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Trans. rate ($\hat\beta$), input cost in LQR design ($\mathbf{R}^{\mathrm{LQR}}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Hospitalized patients ($\hat \mathbf{H}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Time-dependent repropoduction number ($\hat R_t$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} All infected ($\hat \mathbf{L} + \hat \mathbf{P} + \hat \mathbf{I} + \hat \mathbf{A} + \hat \mathbf{H}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Daily new infected ($\hat \beta \, (\hat \mathbf{P} + \hat \mathbf{I} + \hat \delta \hat \mathbf{A}) \, \hat \textbf{S} / \textbf{N}$)' ]
  % [ '\textbf{Plot ' num2str(Cnt) '.} Susceptible people ($\hat \mathbf{S}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} People in the latent phase ($\hat \mathbf{L}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} People in the presympt. phase ($\hat \mathbf{P}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Infected people with sympts. ($\hat \mathbf{I}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Asymptomatic people ($\hat \mathbf{A}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Deceased people ($\hat \mathbf{D}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Recovered and immune people' ]
    % [ '\textbf{Plot ' num2str(Cnt) '.} All recovereds ($\hat \mathbf{R}$)' ]
    % [ '\textbf{Plot ' num2str(Cnt) '.} Immune by vaccination ($\hat \mathbf{U}$)' ]
    };

%%
if args.FigNr == 1
    Cnt(0);

    %% Exp L, P, I, A, H, D
    

    Colors.Color_1 = [0 0.4470 0.7410];
    Colors.Color_2 = [0.8500 0.3250 0.0980];
    Colors.Color_3 = [0.9290 0.6940 0.1250];
    Colors.Color_4 = [0.4940 0.1840 0.5560];
    Colors.Color_5 = [0.4660 0.6740 0.1880];
    Colors.Color_6 = [0.3010 0.7450 0.9330];
    Colors.Color_7 = [0.6350 0.0780 0.1840];

    Perc_Max = 1;
    
    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    Lines = plot(d_recpred,[Exp_L,Exp_P,Exp_I,Exp_A,Exp_H,Exp_D]);
    
    for i = 1:numel(Lines)
        Lines(i).Marker = 'none';
        Lines(i).LineStyle = '-';
        Lines(i).LineWidth = 1.5;
        Lines(i).Color = Colors.(['Color_' num2str(i)]);
    end

    ylim([0,Perc_Max]*r.Par.Np.val/100);
        
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    Leg = legend(Lines,{'$\mathbf{L}$','$\mathbf{P}$','$\mathbf{I}$','$\mathbf{A}$','$\mathbf{H}$','$\mathbf{D}$'},...
        'Interpreter','latex','Location','northwest');
    Leg.NumColumns = 3;
    
    args.ylabel = '';
    args.fname = '_Exp_LPIAHD';
    args.Plotter(r,args,"Show_XLabel",false);
    
    %% beta
    
    % fig = zf_clear_figure(123);
    ax = nexttile;
    hold on
    
    % Std_u_filtered = movmean(Std_u,12);
    Std_u_filtered = Std_u;
    
    Pl = plot(d_recpred(1:end-1),log10(r.R_lqr)/20,'.','LineWidth',1,'Color',[0.8500 0.3250 0.0980]);
    
    Sh = zf_plot_mean_var(d_recpred(1:end-1),Exp_u,Std_u_filtered);
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    % Leg = legend([Sh([2 3 1]) Pl],[ ...
    %     "$1 \times $\,St.dev. $(\approx 68\%)$", ...
    %     "$2 \times $\,St.dev. $(\approx 95\%)$", ...
    %     "Expected value", ...
    %     "$\log_{10}(\mathbf{R}^{\mathrm{LQR}})/20$"
    %     ]);
    Leg = legend([Pl Sh],[ ...
        "$\log_{10}(\mathbf{R}^{\mathrm{LQR}})/20$", ...
        "Exp. value of $\hat \beta$", ...
        "$1 \times $\,St.dev. of $\hat \beta$ $(\approx 68\%)$", ...
        "$2 \times $\,St.dev. of $\hat \beta$ $(\approx 95\%)$"
        ],'Location','northwest');
    Leg.NumColumns = 2;
    
    ylim([0,0.7]);
    
    args.ylabel = '';
    args.fname = '_ExpVar_beta';
    args.Plotter(r,args,"Show_XLabel",false);
    
    %% H
    
    Perc_Max = 0.2;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    Sh = zf_plot_mean_var(d_recpred,Exp_H,Std_H);
    Pl = plot(d_recpred(1:N_rec),H_ref,'Color',[ 0.850 0.325 0.098 ],'Marker','none','LineStyle','-','LineWidth',1.5);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    Leg = legend([Pl Sh(1)],[ ...
        "Official hospitalization ($\mathbf{H}^{\mathrm{Off}}$)", ...
        "Computed expected hospitalization ($\mu^\mathbf{H}$)", ...
        % "$1 \times $\,St.dev. of $\hat \beta$ $(\approx 68\%)$", ...
        % "$2 \times $\,St.dev. of $\hat \beta$ $(\approx 95\%)$"
        ],'Location','northwest');
    Leg.NumColumns = 1;
        
    args.ylabel = '';
    args.fname = '_ExpVar_H';
    args.Plotter(r,args,"Show_XLabel",false);
    

    %% Rt
    
    % fig = zf_clear_figure(123);
    ax = nexttile;
    hold on
            
    Sh = zf_plot_mean_var(d_recpred(1:end-1),r.y_All(3,1:N_rec)',sqrt(r.Sy_diag(3,1:N_rec))');
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    ylim([0,3]);
    
    args.ylabel = '';
    args.fname = '_ExpVar_beta';
    args.Plotter(r,args,"Show_XLabel",false);
    
    %% All infected
    
    Perc_Max = 4;
    
    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    zf_plot_mean_var(r.d_All(1:N_rec),r.y_All(1,1:N_rec)',real(sqrt(r.Sy_diag(1,1:N_rec))'));
    
    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(['~~~~~' Titles{Cnt}],'Interpreter','latex','FontSize',12)
    
    args.ylabel = '';
    args.fname = '_All_inf';
    args.Plotter(r,args,"Show_XLabel",true);
    
    %% Daily new infected
    
    Perc_Max = 0.4;
    
    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    zf_plot_mean_var(r.d_All(1:N_rec),r.y_All(2,1:N_rec)',real(sqrt(r.Sy_diag(2,1:N_rec))'));
    
    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(['~~~~~' Titles{Cnt}],'Interpreter','latex','FontSize',12)
    
    args.ylabel = '';
    args.fname = '_Daily_new';
    args.Plotter(r,args,"Show_XLabel",true);

elseif args.FigNr == 2
    Cnt(6)

    %% S
    
    % ax = nexttile;
    % hold on
    %
    % Sh = zf_plot_mean_var(d_recpred,Exp_S,Std_S);
    % title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    % grid on
    %
    % args.ylabel = '';
    % args.fname = '_ExpVar_L';
    % args.Plotter(r,args,"Show_XLabel",false);
    
    %% L
    
    % fig = zf_clear_figure(1231);
    Perc_Max = 1;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    Sh = zf_plot_mean_var(d_recpred,Exp_L,Std_L);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    args.ylabel = '';
    args.fname = '_ExpVar_L';
    args.Plotter(r,args,"Show_XLabel",false);
    
    %% P
    
    % fig = zf_clear_figure(123);
    Perc_Max = 1;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    zf_plot_mean_var(d_recpred,Exp_P,Std_P);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    args.ylabel = '';
    args.fname = '_ExpVar_P';
    args.Plotter(r,args,"Show_XLabel",false);
    
    %% I
    
    Perc_Max = 1;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    zf_plot_mean_var(d_recpred,Exp_I,Std_I);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    args.ylabel = '';
    args.fname = '_ExpVar_I';
    args.Plotter(r,args,"Show_XLabel",false);
    
    %% A
    
    Perc_Max = 1;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    zf_plot_mean_var(d_recpred,Exp_A,Std_A);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    args.ylabel = '';
    args.fname = '_ExpVar_A';
    args.Plotter(r,args,"Show_XLabel",false);

    %% D
    
    Perc_Max = 0.4;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    zf_plot_mean_var(d_recpred,Exp_D,Std_D);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on
    
    args.ylabel = '';
    args.fname = '_ExpVar_D';
    args.Plotter(r,args,"Show_XLabel",true);
    
    %% R
    
    % fig = zf_clear_figure(123);
    Perc_Max = 70;

    ax = nexttile;
    
    yyaxis right
    ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    ylim([0,Perc_Max])
    
    yyaxis left
    hold on
    Sh = zf_plot_mean_var(d_recpred,Exp_R,Std_R);
    Sh2 = zf_plot_mean_var(d_recpred,Exp_R_All,Std_R_All,[0.8500 0.3250 0.0980]);
    Sh3 = zf_plot_mean_var(d_recpred,Exp_R+Exp_U,sqrt(Std_R.^2 + Std_U.^2),[0.4660 0.6740 0.1880]);

    ylim([0,Perc_Max]*r.Par.Np.val/100);
    
    title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    grid on

    Leg = legend([Sh3(1),Sh2(1),Sh(1)],[ ...
        "Recovered or vaccinated ($\hat\mathbf{R} + \hat\mathbf{U}$)", ...
        "All recovered ($\hat\mathbf{R}^\mathrm{(all)}$)", ...
        "Recovered but not vaccinated ($\hat\mathbf{R}$)"
        ],'Location','northwest');
    Leg.NumColumns = 1;
    % Leg.Box = 'off';
    
    args.ylabel = '';
    args.fname = '_ExpVar_R';
    args.Plotter(r,args,"Show_XLabel",true);
    
    %% U
    % 
    % Perc_Max = 50;
    % 
    % ax = nexttile;
    % 
    % yyaxis right
    % ylabel('\% of population ($\mathbf{N}$)','Interpreter','latex','FontSize',11);
    % ylim([0,Perc_Max])
    % 
    % yyaxis left
    % hold on    
    % zf_plot_mean_var(d_recpred,Exp_U,Std_U);
    % 
    % ylim([0,Perc_Max]*r.Par.Np.val/100);
    % 
    % title(Titles{Cnt},'Interpreter','latex','FontSize',12)
    % grid on
    % 
    % args.ylabel = '';
    % args.fname = '_ExpVar_U';
    % args.Plotter(r,args,"Show_XLabel",false);
    
end

%%

% keyboard
pause(0.5);
fn = append(FilePref,args.fname);

exportgraphics(fig,[ FilePref '_Fig' num2str(args.FigNr) '.png' ])

fname = [ FilePref '_Fig' num2str(args.FigNr) '.pdf' ];
exportgraphics(fig,fname,'ContentType','vector')

fullname = [pwd filesep fname];
fprintf('Grapics exported to:\n%s\n',fullname);
clipboard("copy",fullname);

end

function ret = Cnt(reset)

persistent k

if nargin > 0
    k = reset;
else
    k = k + 1;
end

ret = k;

end


