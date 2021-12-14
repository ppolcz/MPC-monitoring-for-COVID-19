function r = mf_epid_tracking_NMPC(args)
%%
%  File: epid_tracking_NMPC.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study7_tracking
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2021. April 02. (2020b)
%  Minor review on 2021. September 15. (2021a)
%  Reviewed on 2021. October 11. (2021b)
%  Reviewed on 2021. October 12. (2021b)
% 

r.Date_Start_MPC = datetime(2020,08,20);
r.Date_End_REC = datetime(date);
r.Date_End_PRED = r.Date_End_REC + 90;
r.Date_Assume_Beta = r.Date_End_REC + 30;

r.beta_diff_suly = 1e5;
r.beta_ref_suly = 1e3;
r.beta_jovo_suly = 1e4;
r.ref_H_suly = 1e-4;
r.ref_H_rel_suly = true;

r.LastN_beta_ref = 14;
r.beta_jovo_FIX = true;
r.beta_jovo_Ref = 0.65;
r.beta_min = 1/3 * (1 - 0.82);
r.beta_max = 0.7;

r.V_Naponta_jovo = 3000;
r.Oltasok_kozotti_ido = 21;

r.User_prev_sol = false;

r.NL_solver = [];

r = parsepropval('create',r,args);

import casadi.*

%%

% if isempty(r.NL_solver)
%% 

% Index values in the data table (DT)
r.DT_Idx_Start = find(r.DT.Date == r.Date_Start);
r.DT_Idx_Start_MPC = find(r.DT.Date == r.Date_Start_MPC);
r.DT_Idx_End_REC = find(r.DT.Date == r.Date_End_REC);
r.DT_Idx_Last_Available_Data = find(r.DT.Date == r.Date_Last_Available_Data);
r.DT_Idx_Current_Date = daysact(r.Date_Start,datetime(date))+1;

% Check indices
assert(r.DT_Idx_Start == 1);
assert(r.DT_Idx_Last_Available_Data == r.N_Past);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Megj. H_off = H(t*0) [ H(t*1) ... H(t*a.N_rec) ], vagyis a kezdeti      %
% allapot nincs benne.                                                    %
r.H_ref = r.H_ref_All(r.DT_Idx_Start_MPC+1:r.DT_Idx_End_REC);             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dynamical model

s = pf_struct_params(r.Par,"ParamType","exp");
[x,beta,V,f,C,h] = mf_epid_ode_model(r,"ParamType","exp");
[nx,~] = size(x);

dynamics = Function('DT_dynamics',{x,beta,V},{f},{'x','beta','V'},'F');

r.      Fn_dyneq_DET = dynamics;
r.desc__Fn_dyneq_DET = 'Discrete-time deterministic dynamic equations of the system';

%% Load initial condition

if r.Date_Start_MPC > r.DT.Date(1)
    load('rec_and_pred_NMPC_v13_2020-03-01-tol.mat','Rec')
    assert(r.Date_Start == Rec.d(1),['Preliminary reconstruction ' ...
        'should start at the beginnig of the data records.'])
else
    Rec.t = 0;
    Rec.d = r.DT.Date(1);
    Rec.u = 1/3;
    Rec.x = [r.Par.Np.nom-40 , 10 , 10 , 10 , 10 , 0 , 0 , 0 , 0 , 0];
end

r.x0_ref = Rec.x(r.DT_Idx_Start_MPC,:)';
r.u0_ref = Rec.u(r.DT_Idx_Start_MPC);


%% Prediction horizon.

r.N_pref = r.DT_Idx_Start_MPC-1;                    % Preliminarily computed trajectories (past-past)
r.N_rec = r.DT_Idx_End_REC - r.DT_Idx_Start_MPC;    % Reconstruction (past)
r.N_pred = daysact(r.Date_End_REC,r.Date_End_PRED); % Prediction (future)
r.N_recpred = r.N_rec + r.N_pred;

assert(r.N_rec == numel(r.H_ref),'a.N_pred = %d != numel(H_Off) = %d',r.N_rec,numel(r.H_ref));

%% Vaccination model.

r.VV_Available_Data_All = [ r.DT.Vaccinated_1(1:r.N_pref+1+r.N_rec)' zeros(1,r.N_pred) ];
if r.N_pred >= r.Oltasok_kozotti_ido
    r.VV_All = [ zeros(1,r.Oltasok_kozotti_ido) r.DT.Vaccinated_1(1:r.N_pref+1+r.N_rec)' zeros(1,r.N_pred-r.Oltasok_kozotti_ido)+r.V_Naponta_jovo ];
else
    r.VV_All = [ zeros(1,r.Oltasok_kozotti_ido) r.DT.Vaccinated_1(1:r.N_pref+1+r.N_rec+r.N_pred-r.Oltasok_kozotti_ido)' ];
end
r.VV = r.VV_All(r.DT_Idx_Start_MPC:end);

r.VT = table(r.DT.Date(1) + (0:numel(r.VV_All)-1)',r.VV_Available_Data_All',r.VV_All', ...
    'VariableNames',{'Date','V1_Off','V1_Shift_PRED'});

%% Free decision variables along the horizon (column-wise).

helper = Pcz_CasADi_Helper('SX');

xx_fh = @(vars) [ r.x0_ref vars ];
uu_fh = @(vars) [ r.u0_ref vars ];

xx = xx_fh( helper.new_var('x',[nx,r.N_recpred],1,'str','full','lb',0) );
uu = uu_fh( helper.new_var('u',[1,r.N_recpred-1],1,'str','full','lb',r.beta_min,'ub',r.beta_max) );

%% Objective function

w_u = r.beta_ref_suly;
w_du = r.beta_diff_suly;
w_H = r.ref_H_suly; 
w_uf = r.beta_jovo_suly;

if ~r.beta_jovo_FIX
    beta_ref_fh = @(uu) sum(uu(r.N_rec-r.LastN_beta_ref+1:r.N_rec)) / r.LastN_beta_ref;
    r.beta_jovo = beta_ref_fh(uu);    
else
    r.beta_jovo = r.beta_jovo_Ref;
    w_uf = 0;
end
Idx_Assume_Beta = daysact(r.Date_Start_MPC,r.Date_Assume_Beta)+1;
helper.add_obj('Beta_pred',sumsqr(uu(1,Idx_Assume_Beta:end) - r.beta_jovo),w_u)
helper.add_obj('Beta_jovo',(r.beta_jovo_Ref - r.beta_jovo)^2,w_uf*r.N_recpred)

% uu_pred = uu(1,end-a.N_pred:end);
% helper.add_obj('Beta_pred',w_u*sumsqr(uu_pred - uu(a.N_rec)),1)

diff_uu = uu(:,2:end) - uu(:,1:end-1); % <------ minimize l2 norm
helper.add_obj('Beta_rate',sumsqr(diff_uu),w_du)

HH_rec = C*xx(:,2:r.N_rec+1);
if r.ref_H_rel_suly
    % Correction (avoid division by zeros).
    % 2021.11.10. (november 10, szerda), 01:59
    den = r.H_ref; 
    den(den < 1) = min(den(den > 1)); % Avoid division by zero
    H_rel_suly = mean(r.H_ref)./den;
else
    H_rel_suly = 1;
end
helper.add_obj('H_error',(HH_rec - r.H_ref').^2,w_H*H_rel_suly);


%% Equality constraints (initialize)

for k = 1:r.N_recpred
    % fprintf('%d / %d ---> size(uu) = (%d,%d), size(xx) = (%d,%d) \n',k,r.N_recpred,size(uu),size(xx))
    x_kp1 = dynamics(xx(:,k),uu(:,k),r.VV(k));
    helper.add_eq_con( x_kp1 - xx(:,k+1) );
end

r.NL_solver = helper.get_nl_solver;

% end

%%

Solution_found = false;
if isfield(r,'NMPC_sol') && r.User_prev_sol 
    sol_guess = full(r.NMPC_sol.x);
    if ~any(isnan(sol_guess))
        Timer_2xw3 = pcz_dispFunctionName;
        r.NMPC_sol = r.NL_solver.solve([],sol_guess);
        r.NMPC_sol.Elapsed_time = toc(Timer_2xw3);
        r.NMPC_sol.Elapsed_time__desc = 'Previous solution reused';
        pcz_dispFunctionEnd(Timer_2xw3);

        Solution_found = true;
    end
end

if ~Solution_found
    Timer_2xw3 = pcz_dispFunctionName;
    r.NMPC_sol = r.NL_solver.solve();
    r.NMPC_sol.Elapsed_time = toc(Timer_2xw3);
    r.NMPC_sol.Elapsed_time__desc = 'Solution from scratch';
    pcz_dispFunctionEnd(Timer_2xw3);
end

r.NMPC_helper = helper;
[r.NMPC_sol.f,r.NMPC_sol.F,r.NMPC_sol.f_fh] = r.NL_solver.get_obj;

%%

r.xx_val = xx_fh(helper.get_value('x'));
r.uu_val = uu_fh(helper.get_value('u'));
if ~r.beta_jovo_FIX
    r.beta_jovo = beta_ref_fh(r.uu_val);
end

%%

r.t_recpred = 0:r.N_recpred;
r.d_recpred = r.Date_Start_MPC + r.t_recpred';

r.t_All = [ -r.N_pref:-1 , r.t_recpred ];
r.d_All = r.Date_Start_MPC + r.t_All;

r.x_All = [ Rec.x(1:r.DT_Idx_Start_MPC-1,:)' r.xx_val ];
r.u_All = [ Rec.u(1:r.DT_Idx_Start_MPC-1,:)' r.uu_val ];

r.Script_version = sprintf('v15_mf_NMPC_%s',datestr(r.Date_Start_MPC,'yyyy-mm-dd'));

r.y_All = cellfun(@(x,u) {full(h.Fn(x,u))},num2cell(r.x_All,1),num2cell([r.u_All nan]));
r.y_All = horzcat(r.y_All{:});

S_All = r.x_All(1,:)';
L_All = r.x_All(2,:)';
P_All = r.x_All(3,:)';
I_All = r.x_All(4,:)';
A_All = r.x_All(5,:)';
H_All = r.x_All(6,:)';
R_All = r.x_All(7,:)';
dS_All = [0 ; diff(S_All)];
r.Daily_all_Inf = r.y_All(1,:);
r.Daily_new_Inf = r.y_All(2,1:end);
r.Rc_t = r.y_All(3,1:end-1);

if r.Date_Start == r.DT.Date(1)
    Rec.t = r.t_All';
    Rec.d = r.d_All';
    Rec.u = r.u_All';
    Rec.x = r.x_All';
    save('rec_and_pred_NMPC_v15_2020-03-01-tol','Rec')
end

end
