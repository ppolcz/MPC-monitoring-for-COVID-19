function [r] = zf_plotter_opts(r)
%%
%  File: zf_plotter_opts.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%

d_month_All = month(r.d_All);
honap_valtasok_All = [ false , d_month_All(1:end-1) ~= d_month_All(2:end)];

r.Str_Last_Available_Data = sprintf('last available data (%s)',datestr(r.Date_End_REC,'mmm. dd, yyyy'));
r.Str_Current_Date = sprintf('current date (%s)',datestr(datetime(date),'mmm. dd, yyyy'));
r.Str_Date_Start_MPC = sprintf('starting date (%s)',datestr(r.Date_Start_MPC,'mmm. dd, yyyy'));


r.t_Ticks = r.t_All(honap_valtasok_All)';
r.d_Ticks = r.d_All(honap_valtasok_All)';

r.t_Labels = num2cell(datestr(r.d_All(honap_valtasok_All),3),2);

r.t_AxisLabel = sprintf('time [days] $\\in$ [%s ; %s]',...
    datestr(r.Date_Start,'mmm. dd, yyyy')',...
    datestr(r.Date_End_PRED,'mmm. dd, yyyy')');

r.t_Lims = r.t_All([1,end]); % i.e., [-r.N_pref r.N_recpred];
r.d_Lims = r.d_All([1,end]);

% 2021.11.10. (november 10, szerda), 02:02
r.t_Start_MPC = r.t_All(1) + r.DT_Idx_Start_MPC - 1;
r.t_Last_Available_Data = r.t_All(1) + r.DT_Idx_End_REC - 1;
r.t_Current_Date = r.t_All(1) + r.DT_Idx_Current_Date - 1;

r.d_Start_MPC = r.d_All(1) + r.DT_Idx_Start_MPC - 1;
r.d_Last_Available_Data = r.d_All(1) + r.DT_Idx_End_REC - 1;
r.d_Current_Date = r.d_All(1) + r.DT_Idx_Current_Date - 1;

r.Leg_FontSize = 12;
r.Leg_Location = "northoutside";


end