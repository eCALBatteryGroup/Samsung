%% Save Doyle-Fuller-Newman Model Results
%   Created June 13, 2012 by Scott Moura
%   Modified June 26, 2015 by Hector Perez

out.date = date;
out.time = t;
out.cur = I;
out.volt = Volt;
out.soc = SOC;
out.temp = T;

out.c_s_n = c_s_n;
out.c_s_p = c_s_p;
out.c_e = c_e;
out.phi_s_n = phi_s_n;
out.phi_s_p = phi_s_p;
out.i_en = i_en;
out.i_ep = i_ep;
out.phi_e = phi_e;
out.jn = jn;
out.jp = jp;

out.c_ss_n = c_ss_n;
out.c_ss_p = c_ss_p;
out.c_avg_n = c_avg_n;
out.c_avg_p = c_avg_p;
out.SOC = SOC;
out.c_ex = c_ex;
out.eta_n = eta_n;
out.eta_p = eta_p;
out.c_e_0p = c_e_0p;
out.eta_s_Ln = eta_s_Ln;
out.n_Li_s = n_Li_s;
out.n_Li_e = n_Li_e;
out.eta_s_n = eta_s_n;
out.eta_s_p = eta_s_p;   

fileName = input('Save filename? ','s');
save([fileName '.mat'],'out');