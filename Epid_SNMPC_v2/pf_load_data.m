function r = pf_load_data(r)
arguments
    r = struct
end
%%
%  File: pf_load_data.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study7_tracking
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2021. March 31. (2020b)
%  Major review on 2021. September 20. (2021a)
%  Reviewed on 2021. October 12. (2021b)
%

%% Load xls

xls_name = 'data_All_2019-03-31.xls';
Wiki_angol = readtable(xls_name, 'ReadVariableNames', true, 'Sheet', 'Wiki_angol');

%%

Rt_data = readtable(xls_name, 'ReadVariableNames', true, 'Sheet', 'Rt');

Rt_Good = ~isnan(Rt_data.Rt);
N_Rt = numel(Rt_Good);

t = 0:N_Rt-1;
t_Good = t(Rt_Good);
Rt_Good = Rt_data.Rt(Rt_Good);

t_Relevand = t_Good(1):t_Good(end);
Rt_Relevant = interp1(t_Good,Rt_Good,t_Relevand,'linear');
Rt_ma7 = movmean(Rt_Relevant,7);

r.AtloTeam.t_Rt = t_Relevand;
r.AtloTeam.d_Rt = Rt_data.Date(1) + t_Relevand;
r.AtloTeam.Rt = Rt_Relevant;
r.AtloTeam.Rt_ma7 = Rt_ma7;

%% Date variables

Start_Date = datetime(2020,03,01);

Start_Date_of_Data = datetime(Wiki_angol.Date{1},'InputFormat','yyyy.MM.dd.');
End_Date_of_Data = datetime(Wiki_angol.Date{end},'InputFormat','yyyy.MM.dd.');

%% Simulation time span

N = days(End_Date_of_Data - Start_Date);
t = (0:N)';
t_date = Start_Date + t;

%{
    table(t,t_date)
%}

%% Prepare hospitalization and vaccination time series

Padding = zeros(days(Start_Date_of_Data - Start_Date),1);

H = [
    Padding
    Wiki_angol.Hospitalized
    ];

V_first = [
    Padding
    Wiki_angol.ATLO_Elso
    ];

V_second = [
    Padding
    Wiki_angol.ATLO_Masodik
    ];

H(isnan(H) | H < 0) = 0;
V_first(isnan(V_first)) = 0;
V_second(isnan(V_second)) = 0;


%% Hogy jobban lassunk
% 2021.06.04. (június  4, péntek), 10:05
% 2021.09.20. (szeptember 20, hétfő), 14:23

r.DT = table(t_date,t,H(1:N+1),V_first(1:N+1),V_second(1:N+1),'VariableNames',{'Date' 'Day' 'Hospitalized' 'Vaccinated_1' 'Vaccinated_2'});

r.N_Past = size(r.DT,1);

% Filtered data series of hospitalized patients.
r.H_ref_All = movmean(r.DT.Hospitalized,7);

%% Dates

r.Date_Start = r.DT.Date(1);
r.Date_Last_Available_Data = r.DT.Date(end);
r.Date_End_REC = r.Date_Last_Available_Data;

end
