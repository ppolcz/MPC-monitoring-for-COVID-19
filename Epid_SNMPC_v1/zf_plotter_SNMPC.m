function zf_plotter_SNMPC(r,args)
arguments
    r
    args.Layout = {5 2};
    args.FigSize = [915 1182];
    args.FigNr = 8081;
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

FilePref = sprintf(['results/SNMPC%s/%s_Pred-%s'],...
    args.Version,datestr(date,29),datestr(r.Date_End_REC,29));
mkdir(fileparts(FilePref))

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
  % [ '\textbf{Plot ' num2str(Cnt) '.} Susceptible people ($\hat \mathbf{S}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} People in the latent phase ($\hat \mathbf{L}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} People in the presympt. phase ($\hat \mathbf{P}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Infected people with sympts. ($\hat \mathbf{I}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Asymptomatic people ($\hat \mathbf{A}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Hospitalized patients ($\hat \mathbf{H}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Recovered people ($\hat \mathbf{R}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Deciesed people ($\hat \mathbf{D}$)' ]
    [ '\textbf{Plot ' num2str(Cnt) '.} Immune by vaccination ($\hat \mathbf{U}$)' ]
    };
Cnt(0);

%% Exp L, P, I, A, H, D

ax = nexttile;
hold on
Lines = plot(d_recpred,[Exp_L,Exp_P,Exp_I,Exp_A,Exp_H,Exp_D]);

% Lines = [];
% for i = [2:6 8]
%     Lines = [Lines zf_plot_mean_var(sim_t,Mxx_val(i,:)',sqrt(Sxx_diag(i,:)'),pcz_get_plot_colors([],i-1))];
% end

title(Titles{Cnt},'Interpreter','latex','FontSize',12)
Leg = legend(Lines,{'$\mathbf{L}$','$\mathbf{P}$','$\mathbf{I}$','$\mathbf{A}$','$\mathbf{H}$','$\mathbf{D}$'},...
    'Interpreter','latex','Location','northwest');
Leg.NumColumns = 3;

args.ylabel = '';
args.fname = '_Exp_LPIAHD';
args.Plotter(r,args,"Show_XLabel",false);

%% beta

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
    "Expected value of $\hat \beta$", ...
    "$1 \times $\,St.dev. of $\hat \beta$ $(\approx 68\%)$", ...
    "$2 \times $\,St.dev. of $\hat \beta$ $(\approx 95\%)$"
    ]);
Leg.NumColumns = 2;

ylim([0,0.6]);

args.ylabel = '';
args.fname = '_ExpVar_beta';
args.Plotter(r,args,"Show_XLabel",false);

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

ax = nexttile;
hold on

Sh = zf_plot_mean_var(d_recpred,Exp_L,Std_L);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_L';
args.Plotter(r,args,"Show_XLabel",false);

%% P

ax = nexttile;
hold on

zf_plot_mean_var(d_recpred,Exp_P,Std_P);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_P';
args.Plotter(r,args,"Show_XLabel",false);

%% I

ax = nexttile;
hold on

zf_plot_mean_var(d_recpred,Exp_I,Std_I);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_I';
args.Plotter(r,args,"Show_XLabel",false);

%% A

ax = nexttile;
hold on

zf_plot_mean_var(d_recpred,Exp_A,Std_A);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_A';
args.Plotter(r,args,"Show_XLabel",false);

%% H

ax = nexttile;
hold on

Line1 = zf_plot_mean_var(d_recpred,Exp_H,Std_H);
Line2 = plot(d_recpred(1:N_rec),H_ref,'Color',[ 0.850 0.325 0.098 ]);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

legend([Line1(1) Line2], { ...
    '($\mu^\mathbf{H}$) Computed hospitatization', ...
    '($\mathbf{H}^{\mathrm{Off}}$) Official data' ...
    },'Location','northwest')

r.Leg_Location = "best";

args.fname = '_ExpVar_H';
args.Plotter(r,args,"Show_XLabel",false);

%% R

ax = nexttile;
hold on

zf_plot_mean_var(d_recpred,Exp_R,Std_R);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_R';
args.Plotter(r,args,"Show_XLabel",false);

%% D

ax = nexttile;
hold on

zf_plot_mean_var(d_recpred,Exp_D,Std_D);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_D';
args.Plotter(r,args);

%% U

ax = nexttile;
hold on

zf_plot_mean_var(d_recpred,Exp_U,Std_U);
title(Titles{Cnt},'Interpreter','latex','FontSize',12)
grid on

args.ylabel = '';
args.fname = '_ExpVar_U';
args.Plotter(r,args);


%%

%{

ax = nexttile;

N_prefrec = r.N_pref + r.N_rec;

hold on
plot(d_All,r.x_All([2:6,8],:))
plot(d_All(1:N_prefrec),r.H_ref_All(1:N_prefrec)','r','LineWidth',1.5)
% plot(sim_t(1:end-1),u_val_pred*1e5,'k','LineWidth',1.5)

YMax = 90000;
ax.YLim = [0 YMax];

title('($\mathbf{L}$--$\mathbf{H},\mathbf{D}$) Estimated trajectories for fixed parameters',...
    'Interpreter','latex','FontSize',11)
Leg = legend('$\mathbf{L}$','$\mathbf{P}$','$\mathbf{I}$','$\mathbf{A}$','$\mathbf{H}$','$\mathbf{D}$',...
    '$\mathbf{H}^{\mathrm{Off}}$ (Official hospitalization data)',...
    ...'($\beta \times 10^5$) Transmission rate',...
    'Interpreter','latex','Location','northwest');

Leg.NumColumns = 2;

yticks((0:13)*10000)

args.ylabel = '';
args.fname = '_Fig1_LPIAHD';
args.Plotter(r,args);

%}

%%

% keyboard
pause(0.5);
fn = append(FilePref,args.fname);

exportgraphics(fig,[ FilePref '_Fig' num2str(args.FigNr) '.png' ])

% for i = 1:numel(Ax)
%     exportgraphics(Ax(i),[ FilePref '_Ax' num2str(i) '.png' ])
% end

fname = [ FilePref '_Fig' num2str(args.FigNr) '.pdf' ];
exportgraphics(fig,fname,'ContentType','vector')

fullname = [pwd filesep fname];
fprintf('Grapics exported to:\n%s\n',fullname);
clipboard("copy",fullname);

end

function ret = Cnt(reset)

persistent k

if nargin > 0
    k = 0;
else
    k = k + 1;
end

ret = k;

end


