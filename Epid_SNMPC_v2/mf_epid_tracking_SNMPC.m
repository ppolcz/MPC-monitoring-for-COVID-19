function r = mf_epid_tracking_SNMPC(r)
%%
%  File: mf_epid_tracking_SNMPC.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2021. November 17. (2021b)
%

import casadi.*

He = @(A) A + A';

params = r.params;
Sp_val = r.Sp;
Mp_val = r.Mp;

N_rec = r.N_rec;
N_recpred = r.N_recpred;

H_ref = r.H_ref;
VV = r.VV;

Mxx_guess = r.Mxx_val;
Muu_guess = r.Muu_val;
Sxx_guess = r.Sxx_val;
Sxxp_guess = r.Sxxp_val;
KK = r.KK_val;

%%%
% Dynamical model

[x,beta,V,f,C] = mf_epid_ode_model(r,"ParamType","sym");
[nx,~] = size(x);
np = numel(params);

% Presumed variance matrix of the initial state.
Sx0 = blkdiag(7,eye(7),0.01);

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
    
    full(f_Ak([1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    full(f_Bk([1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    full(f_Ek([1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    full(f_ck([1;1;1;1;1;1;1;1],1,[1;1;1;1;1;1;1;1;1;1],1))
    
%}

% Symbolic placeholder variables: mu_u(k), mu_x(k), Sigma_x(k), Sigma_xp(k)
% Mu = SX.sym('mu_u_free',size(beta));
% Mx = SX.sym('mu_x_free',size(x));
[Sx,Sx_half,Fn_Sx_vec] = Pcz_CasADi_Helper.create_sym('Sigma_x',nx,'str','sym','r3','f_vec');
[Sxp,Fn_Sxp_vec] = Pcz_CasADi_Helper.create_sym('Sigma_xp',[nx np],'str','full','r2','f_vec');

% Symbolic placeholder variables: mu_u(k+1), Sigma_xp(k+1)
beta_pp = SX.sym('beta_pp',size(beta));
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

Exp_du = (beta_pp-beta)^2 + K_pp*Sx_pp*K_pp' + K*Sx*K' - He( K_pp*( (A-B*K)*Sx + E*Sxp' )*K' );
Fn_Exp_du = Function('Fn_Exp_du', ...
    {x,beta,beta_pp,V,params,Sxp,Sx_half,Sx_pp_half,K,K_pp}, ...
    {Exp_du}, ...
    {'Mx','Mu','Mu_pp','V','Mp','Sxp','Sx','Sx_pp','K','K_pp'}, ...
    {'Exp_du'});

%% Free decision variables along the horizon (column-wise).

helper = Pcz_CasADi_Helper('SX');

Mxx_fh = @(vars) [ Mxx_guess(:,1) vars ];
Muu_fh = @(vars) [ Muu_guess(:,1) vars ];
Sxx_fh = @(vars) [ Sxx_guess(1) vars ];
Sxxp_fh = @(vars) [ Sxxp_guess(1) vars ];

Mxx = Mxx_fh( helper.new_var('mu_x',[nx,N_recpred],1,'str','full','lb',0) );
Muu = Muu_fh( helper.new_var('mu_u',[1,N_recpred-1],1,'str','full','lb',1/3*(1 - 0.82),'ub',1) );
Sxx = Sxx_fh( helper.new_var('Sigma_x',nx,N_recpred,'str','sym') );
Sxxp = Sxxp_fh( helper.new_var('Sigma_xp',[nx np],N_recpred,'str','full') );

%% Equality constraints (initialize)

status = PStatus(N_recpred,'Construct prediction model');
for k = 1:N_recpred
    
    % f_mu_x_kp1:(x_hat[8],u_hat,V,mu_p[3],x[8],u)->(mu_x_kp1[8])
    mu_x_kp1 = Fn_Mx_pp(Mxx(:,k),Muu(:,k),VV(k),Mp_val);
    
    % f_Sigma_xp_kp1:(x_hat[8],u_hat,V,mu_p[3],Sigma_xp[8x3])->(Sigma_xp_kp1[8x3])
    Sigma_xp_kp1 = Fn_Sxp_pp(Mxx(:,k),Muu(:,k),VV(k),Mp_val,Sxxp{k},KK(k,:));
    
    % f_Sigma_x_kp1:(x_hat[8],u_hat,V,mu_p[3],Sigma_xp[8x3],Sigma_x[8x8,36nz])->(Sigma_x_kp1[8x8])
    Sigma_x_kp1 = Fn_Sx_pp(Mxx(:,k),Muu(:,k),VV(k),Mp_val,Sxxp{k},Sxx{k},KK(k,:));

    helper.add_eq_con( mu_x_kp1 - Mxx(:,k+1) );
    helper.add_eq_con( Fn_Sxp_vec( Sigma_xp_kp1 - Sxxp{k+1} ) );
    helper.add_eq_con( Fn_Sx_vec( Sigma_x_kp1 - Sxx{k+1} ) );

    status.progress(k);
end

%%  Objective function

% Weights:
w_du = r.beta_diff_suly;
w_H = r.ref_H_suly; 

HH_rec = C*Mxx(:,2:N_rec+1);
if r.ref_H_rel_suly
    % Correction (avoid division by zeros).
    % 2021.11.10. (november 10, szerda), 01:59
    den = r.H_ref; 
    den(den < 1) = min(den(den > 1)); % Avoid division by zero
    H_rel_suly = mean(r.H_ref)./den;
else
    H_rel_suly = ones(size(r.H_ref));
end

status = PStatus(N_recpred-1,'Objective function''s `du` term');
for k = 1:N_recpred-1
    E_du_k = Fn_Exp_du(Mxx(:,k),Muu(:,k),Muu(:,k+1),VV(k),Mp_val,Sxxp{k},Sxx{k},Sxx{k+1},KK(k,:),KK(k+1,:));
    helper.add_obj('Beta_rate',E_du_k,w_du);
    status.progress(k);
end

status = PStatus(N_rec+1,'Objective function''s `H - H_ref` term');
helper.add_obj('H_error',(HH_rec - H_ref').^2,w_H*H_rel_suly)
status.progress(1);
for k = 1:N_rec
    helper.add_obj('H_error',C*Sxx{k+1}*C',w_H*H_rel_suly(k));
    status.progress(k+1);
end

% for k = 1:N_recpred
%     helper.add_obj('x_error',w_u*(Muu(k)-Muu_guess(k))^2 + w_u*KK(k,:)*Sxx{k}*KK(k,:)',1);
% end

% Objective function's `predictive u` term:
% for k = Idx_Assume_Beta:N_recpred
%     helper.add_obj('Beta_pred',w_u*(Muu(k)-beta_Jovo)^2 + w_u*KK(k,:)*Sxx{k}*KK(k,:)',1)
% end

% Get solver object
solver = helper.get_nl_solver;

Fn_var = helper.gen_var_mfun;
r.SNMPC_guess.x = full(Fn_var(Mxx_guess(:,2:end),Muu_guess(:,2:end),Sxx_guess{2:end},Sxxp_guess{2:end}));
[r.SNMPC_guess.f,r.SNMPC_guess.F] = helper.get_obj(r.SNMPC_guess.x);

%% Solve the QP problem

r.SNMPC_helper = helper;

Timer_2xw3 = pcz_dispFunctionName;
r.SNMPC_sol = solver.solve([],r.SNMPC_guess.x);
r.SNMPC_sol.Elapsed_time = toc(Timer_2xw3);
pcz_dispFunctionEnd(Timer_2xw3);

%%

r.Muu_val = Muu_fh(helper.get_value('mu_u'));
r.Mxx_val = Mxx_fh(helper.get_value('mu_x'));
r.Sxx_val = Sxx_fh(helper.get_value('Sigma_x'));

Sxx_diag_cell = cellfun(@(Sigma) {diag(Sigma)},r.Sxx_val);
r.Sxx_diag = [ Sxx_diag_cell{:} ];

Suu_val = zeros(1,N_recpred);
for k = 1:N_recpred
    Suu_val(k) = KK(k,:) * r.Sxx_val{k} * KK(k,:)';
end

r.Suu_val = Suu_val;
r.KK_val = KK;

end

%%

% epid_plot_SNMPC_LTV_with_fdbk
