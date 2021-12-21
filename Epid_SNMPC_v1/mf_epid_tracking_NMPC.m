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

r.w_diff_beta = 1e4;
r.w_ref_H = 1e-4;
r.w_ref_H_is_relative = true;

r.beta_min = 1/3 * (1 - 0.82);
r.beta_max = 0.7;

% Average number of days, after which the first dose makes effect
r.Tv = 21;

r.NL_solver = [];

r = parsepropval('create',r,args);

import casadi.*

%%

% Index values in the data table (DT)
r.DT_Idx_Start_MPC = find(r.DT.Date == r.Date_Start_MPC);
r.DT_Idx_End_REC = find(r.DT.Date == r.Date_End_REC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Official hospitalization data:                                          %
% H_off = H(t*0) [ H(t*1) ... H(t*a.N_rec) ],                             %
% (The initial value H(t*0) is not considered.)                           %
r.H_ref = r.H_ref_All(r.DT_Idx_Start_MPC+1:r.DT_Idx_End_REC);             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dynamical model

s = pf_struct_params(r.Par,"ParamType","exp");
[x,beta,V,f,C,h] = mf_epid_ode_model(r,"ParamType","exp");
[nx,~] = size(x);

dynamics = Function('DT_dynamics',{x,beta,V},{f},{'x','beta','V'},'F');

r.      Fn_dyneq_DET = dynamics;
r.desc__Fn_dyneq_DET = 'Discrete-time deterministic dynamic equations of the system';

%% Initial condition

r.x0_ref = [s.Np-40 , 10 , 10 , 10 , 10 , 0 , 0 , 0 , 0 , 0]';
r.u0_ref = 1/3;

%% Prediction horizon.

% Reconstruction (past):
r.N_rec = r.DT_Idx_End_REC - r.DT_Idx_Start_MPC; 

assert(r.N_rec == numel(r.H_ref),'a.N_pred = %d != numel(H_Off) = %d',r.N_rec,numel(r.H_ref));

%% Vaccination model.

r.VV_Available_Data = r.DT.Vaccinated_1(1:r.N_rec+1)';
r.VV = [ zeros(1,r.Tv) r.DT.Vaccinated_1(1:r.N_rec+1-r.Tv)' ];
r.VT = table(r.DT.Date(1) + (0:r.N_rec)',r.VV_Available_Data',r.VV', ...
    'VariableNames',{'Date','V1_Off','V1_Shifted'});

%% Free decision variables along the horizon (column-wise).

helper = Pcz_CasADi_Helper('SX');

xx_fh = @(vars) [ r.x0_ref vars ];
uu_fh = @(vars) [ r.u0_ref vars ];

xx = xx_fh( helper.new_var('x',[nx,r.N_rec],1,'str','full','lb',0) );
uu = uu_fh( helper.new_var('u',[1,r.N_rec-1],1,'str','full','lb',r.beta_min,'ub',r.beta_max) );

%% Objective function

diff_uu = uu(:,2:end) - uu(:,1:end-1); % <------ minimize l2 norm
helper.add_obj('Beta_rate',sumsqr(diff_uu),r.w_diff_beta)

HH_rec = C*xx(:,2:r.N_rec+1);
if r.w_ref_H_is_relative
    % Correction (avoid division by zeros).
    % 2021.11.10. (november 10, szerda), 01:59
    den = r.H_ref; 
    den(den < 1) = min(den(den > 1)); % Avoid division by zero
    w_ref_H_relative = mean(r.H_ref)./den;
else
    w_ref_H_relative = 1;
end
helper.add_obj('H_error',(HH_rec - r.H_ref').^2,r.w_ref_H*w_ref_H_relative);


%% Equality constraints (initialize)

for k = 1:r.N_rec
    x_kp1 = dynamics(xx(:,k),uu(:,k),r.VV(k));
    helper.add_eq_con( x_kp1 - xx(:,k+1) );
end

r.NL_solver = helper.get_nl_solver;

%%

Timer_2xw3 = pcz_dispFunctionName;
r.NMPC_sol = r.NL_solver.solve();
r.NMPC_sol.      Elapsed_time = toc(Timer_2xw3);
r.NMPC_sol.desc__Elapsed_time = 'Solution from scratch';
pcz_dispFunctionEnd(Timer_2xw3);

r.NMPC_helper = helper;
[r.NMPC_sol.f,r.NMPC_sol.F,r.NMPC_sol.f_fh] = r.NL_solver.get_obj;

%%

r.t_val = 0:r.N_rec;
r.d_val = r.Date_Start_MPC + r.t_val;
r.x_val = xx_fh(helper.get_value('x'));
r.u_val = uu_fh(helper.get_value('u'));

r.y_val = cellfun(@(x,u) {full(h.Fn(x,u))},num2cell(r.x_val,1),num2cell([r.u_val nan]));
r.y_val = horzcat(r.y_val{:});

r.Daily_all_Inf = r.y_val(1,:);
r.Daily_new_Inf = r.y_val(2,1:end);
r.Rc_t = r.y_val(3,1:end-1);

%%

r.Script_version = sprintf('v1_mf_NMPC_%s',datestr(r.Date_Start_MPC,'yyyy-mm-dd'));

end
