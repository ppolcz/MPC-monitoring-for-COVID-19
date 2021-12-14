%%
%  File: Main0_README.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study13_SNMPC_LTV_delta
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. October 12. (2021b)
%
%  Az időskála három NAGY részre van osztva: 1. a tartomány, amit korábban
%  rekonstruáltam és már nem akarom változtatni; 2. a tartomány, amit
%  múltnak tekintek, ahol a kórházban ápoltak számát rekonstruálni
%  szeretném; 3. a tartomány, amit jövőnek tekintek és azt előre szeretném
%  jelezni.
% 
%   (Ha Date_End_REC == Date_Last_Available_Data, akkor a JOVO == PRED.)  ┊
%       ┕━━━━━━━━━━┥    ┕━━━━━━━━━━━━━━━━━━━━━━┵─────╮                    ┊
%                  ╰────────────────────────────╮    ┊                    ┊
% Date_Start (01-Mar-2020)                      ┊    ┊                    ┊    
% ┕━━━━━━━━┵─╮ Date_Start_MPC (eg. 20-Aug-2020) ┊    ┊                    ┊
%            ┊ ┕━━━━━━━━━━━━┵─╮                 ┊    ┊                    ┊
%            ┊                ┊                 ┊    ┊                    ┊
%            ┊<---- Official_data_table (MULT) -┊--->┊<----- JOVO ------->┊
%            ┊                ┊                 ┊    ┊                    ┊
%            ┊<---- PREF -----┊-->┊<--- REC --->┊<----- PRED ------------>┊
%            ┊                ┊   ┊             ┊    ┊                    ┊
%   t_All  = [ -N_pref ... -1 ┊ 0 ┊ 1 ... N_rec ┊ N_rec+1 ..... N_recpred ]
%   t_recpred =               [ 0 ┊             ┊    ┊          N_recpred ]
%            ┊                ┊   ┊             ┊    ┊                    ┊
%   DT_Idx:  [  1  .......... ┊ ~ ┊ ......... ~ ┊  ~ ]                    ┊
%       DT_Idx_Start_MPC────────╯ ┊           │ ┊  │ ┊                    ┊
%       DT_Idx_End_REC────────────────────────╯ ┊  │ ┊                    ┊
%       DT_Idx_Last_Available_Data == N_Past───────╯ ┊                    ┊
%                             ┊   ┊             ┊    ┊                    ┊
%   H_ref_All:<---------------┊---┊-------------┊--->┊  (filtered)        ┊
%   H_ref (Recon. with MPC):  ┊   ┊<----------->┊    ┊  (filtered)        ┊
%  ┍━━━━━━━━━━━━━━━━━━━┑      ┊   ┊             ┊    ┊                    ┊
%  ┆ Time span of MPC: ┆      ┝━━━┿━━━━━━━━━━━━━┿━━━━┿━━━━━━━━━━━━━━━━━━━━┥
%  ┕━━━━━━━━━━━━━━━━━━━┙      ┊   ┊             ┊    ┊                    ┊
%   xx:                       [x0 ┊ CasADi symvars 1:N_recpred      ~   ~ ]
%   uu:                       [u0 ┊ CasADi symvars 1:N_recpred-1    ~ ]   ┊
%   Prescr. smoothness for u: ┊<--┊-------------┊----┊-----┊--------->┊   ┊
%   Prescr. reference for u:  ┊   ┊                        ┊<-------->┊   ┊
%   Indexing of uu and xx     [ 1 ┊  Date_Assume_Beta      ┊ ~      ~ ] ~ ┊
%                                    ┕━━━━━━━━━━━━━━┵──────╯ │      │   │ ┊
%                                     Idx_Assume_Beta────────╯      │   │ ┊
%                                              N_recpred────────────╯+1─╯ ┊

