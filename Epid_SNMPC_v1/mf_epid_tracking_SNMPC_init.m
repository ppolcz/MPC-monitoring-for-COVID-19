function r = mf_epid_tracking_SNMPC_init(r)
%%
%  File: mf_epid_tracking_SNMPC_init.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2021. November 16. (2021b)
% 
% Forked from:
% 
%  File: epid_tracking_SNMPC_LTV.m
%  Author: Peter Polcz (ppolcz@gmail.com)
%  Reviewed: 2021. October 12. (2021b)
%

default = struct;
r = parsepropval('create',default,r);

%%

import casadi.*

He = @(A) A + A';

%%%
% Collect data from the NMPC solution

% Symbolic parameter objects 
params = r.params;

% Mean and variance of the parameters
Mp_val = r.Mp;
Sp_val = r.Sp;

N_rec = r.N_rec;
N_recpred = r.N_recpred;

H_ref = r.H_ref;
VV = r.VV;

xx_val_NMPC = r.xx_val;
uu_val_NMPC = r.uu_val;

%%%
% Dynamical model

[x,beta,V,f,C,h,hvar] = mf_epid_ode_model(r,"ParamType","sym");
[nx,~] = size(x);
np = numel(params);

% Presumed variance matrix of the initial state.
Sx0 = blkdiag(7,eye(7),0.01,1);

A = jacobian(f,x);
B = jacobian(f,beta);
E = jacobian(f,params);
c = f - B*beta - A*x - E*params;

Fn_A = Function('Fn_A',{x,beta,params,V},{A},...
    {'x','u','p','V'},'A');
Fn_B = Function('Fn_B',{x,beta,params,V},{B},...
    {'x','u','p','V'},'Bk');
Fn_E = Function('Fn_E',{x,beta,params,V},{E},...
    {'x','u','p','V'},'E');
Fn_c = Function('Fn_c',{x,beta,params,V},{c},...
    {'x','u','p','V'},'c');

%{
    GYORS TEST:
    
    full(Fn_A([1;1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    full(Fn_B([1;1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    full(Fn_E([1;1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    full(Fn_c([1;1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    
%}

% Symbolic placeholder variables: mu_u(k), mu_x(k), Sigma_x(k), Sigma_xp(k)
Mu = SX.sym('mu_u_free',size(beta));
Mx = SX.sym('mu_x_free',size(x));
[Sx,Sx_half,Fn_Sx_vec] = Pcz_CasADi_Helper.create_sym('Sigma_x',nx,'str','sym','r3','f_vec');
[Sxp,Fn_Sxp_vec] = Pcz_CasADi_Helper.create_sym('Sigma_xp',[nx np],'str','full','r2','f_vec');

% Symbolic placeholder variables: mu_u(k+1), Sigma_xp(k+1)
Mu_pp = SX.sym('mu_u_kp1',size(beta));
[Sx_pp,Sx_pp_half] = Pcz_CasADi_Helper.create_sym('Sigma_x_pp',nx,'str','sym');

% Symbolic placeholder variables for K(k) and K(k+1)
K = SX.sym('K_k',size(B'));
K_pp = SX.sym('K_kp1',size(B'));

Fn_Mx_pp = Function('Fn_Mx_pp',...
    {x,beta,V,params},...
    {f},...
    {'Mx','Mu','v','Mp'},...
    {'Mx_pp'});

Fn_Sxp_pp = Function('Fn_Sxp_pp',...
    {x,beta,V,params,Sxp,K},...
    {(A-B*K)*Sxp + E*Sp_val},...
    {'Mx','Mu','v','Mp','Sxp','K'},...
    {'Sxp_pp'});

Fn_Sx_pp = Function('Fn_Sx_pp',...
    {x,beta,V,params,Sxp,Sx_half,K},...
    {(A-B*K)*Sx*(A-B*K)' + E*Sp_val*E' + He( (A-B*K)*Sxp*E' )},...
    {'Mx','Mu','v','Mp','Sxp','Sx','K'},...
    {'Sx_pp'});

Exp_du = (Mu_pp-Mu)^2 + K_pp*Sx_pp*K_pp' + K*Sx*K' - He( K_pp*( (A-B*K)*Sx + E*Sxp' )*K' );
Fn_Exp_du = Function('Fn_Exp_du', ...
    {x,beta,V,params,Sxp,Sx_half,Sx_pp_half,K,K_pp}, ...
    {Exp_du}, ...
    {'Mx','Mu','V','Mp','Sxp','Sx','Sx_pp','K','K_pp'}, ...
    {'Exp_du'});


r.      Fn_dyneq_STO = Fn_Mx_pp;
r.desc__Fn_dyneq_STO = 'Discrete-time dynamic equations of the expected state';

%%%
% Compute variances

idx = 1:6;

Mxx_pp_error = zeros(N_recpred,1);
KK = zeros(N_recpred,nx);
Sxx = [ {Sx0} cell(1,N_recpred) ];
Sxxp = [ {zeros(nx,np)} cell(1,N_recpred) ];
Suu = zeros(1,N_recpred);
R_lqr = zeros(1,N_recpred);

% 2021.12.14. (december 14, kedd), 16:08
Sxup = cell(1,N_recpred);
Sxup_fh = @(Sx,K,Sxp) ...
    blkdiag([eye(nx) -K'],eye(np))' * [ 
        Sx   Sxp
        Sxp' Sp_val
    ] * blkdiag([eye(nx) -K'],eye(np));
% Sxup_fh = @(Sx,K,Sxp) [
%      Sx   -Sx*K'    Sxp
%     -K*Sx  K*Sx*K' -K*Sxp
%     Sxp'  -Sxp*K    Sp_val
%     ];

Nr_MaxIt = 100;

status = PStatus(N_recpred,'Compute variances for initial guess');
for k = 1:N_recpred
    Mx_val = xx_val_NMPC(:,k);
    Mu_val = uu_val_NMPC(k);

    % Prediction error.
    % f_mu_x_kp1:(x_hat[8],u_hat,V,mu_p[3],x[8],u)->(mu_x_kp1[8])
    mu_x_kp1 = Fn_Mx_pp(Mx_val,Mu_val,VV(k),Mp_val);
    Mxx_pp_error(k) = norm(full(mu_x_kp1) - xx_val_NMPC(:,k+1))^2;
    
    % Compute feedback gain
    A_num = full(Fn_A(Mx_val,Mu_val,Mp_val,VV(k)));
    B_num = full(Fn_B(Mx_val,Mu_val,Mp_val,VV(k)));    
    E_num = full(Fn_E(Mx_val,Mu_val,Mp_val,VV(k)));

    Std_u = 100;
    It = 1;
    R_lqr_fh = @(i) 2^(i-1);
    while It < Nr_MaxIt && (Mu_val-2*Std_u <= 0 || 1 <= Mu_val+2*Std_u)
        KK(k,idx) = dlqr(A_num(idx,idx),B_num(idx,:),eye(numel(idx)),R_lqr_fh(It));    
        Suu(k) = KK(k,:) * Sxx{k} * KK(k,:)';
        Std_u = sqrt(Suu(k));
        It = It + 1;
    end
    R_lqr(k) = R_lqr_fh(It-1);

    % Joint variance of the actual state, input, and parameter
    % 2021.12.14. (december 14, kedd), 16:09
    Sxup{k} = Sxup_fh(Sxx{k},KK(k,:),Sxxp{k});
    
    % Compute state-parameter covariance
    % f_Sigma_xp_kp1:(x_hat[8],u_hat,V,mu_p[3],Sigma_xp[8x3])->(Sigma_xp_kp1[8x3])
    Sxxp{k+1} = full(Fn_Sxp_pp(Mx_val,Mu_val,VV(k),Mp_val,Sxxp{k},KK(k,:)));
    
    % Compute state variance
    % f_Sigma_x_kp1:(x_hat[8],u_hat,V,mu_p[3],Sigma_xp[8x3],Sigma_x[8x8,36nz])->(Sigma_x_kp1[8x8])
    Sxx{k+1} = full(Fn_Sx_pp(Mx_val,Mu_val,VV(k),Mp_val,Sxxp{k},Sxx{k},KK(k,:)));

    status.progress(k);
end

Sxup = [Sxup {Sxup_fh(Sxxp{end},0*KK(end,:),Sxxp{end})}];

r.Sy_diag = cellfun(@(x,u,S) {full(hvar.Fn(x,u,Mp_val,S))}, ...
    num2cell(r.x_All,1), ...
    num2cell([r.u_All nan]), Sxup);
r.Sy_diag = horzcat(r.Sy_diag{:});


r.Muu_val = uu_val_NMPC;
r.Mxx_val = xx_val_NMPC;
r.Sxx_val = Sxx;
r.Sxxp_val = Sxxp;
r.Sxup_val = Sxup;
r.Mx_pp_error = Mxx_pp_error;
r.Suu_val = Suu;
r.KK_val = KK;
r.R_lqr = R_lqr;

Sxx_diag_cell = cellfun(@(Sigma) {diag(Sigma)},r.Sxx_val);
r.Sxx_diag = [ Sxx_diag_cell{:} ];

end
