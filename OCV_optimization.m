clc
clear all
close all
params_NMC_Samsung
load Vmean2_Co20_end2p6.mat

cell_capacity = 2.08;%Ah
max_Capacity_P = 4.34;%Ah
max_Capacity_N = 3.63;%Ah

max_concentration_P = (max_Capacity_P*3600)/(p.epsilon_s_p*p.L_p*p.Area*p.Faraday);%mol/m3
max_concentration_N = (max_Capacity_N*3600)/(p.epsilon_s_n*p.L_n*p.Area*p.Faraday);%mol/m3

P_stoichiometry_diff_0_100_SOC = cell_capacity/max_Capacity_P
N_stoichiometry_diff_0_100_SOC = cell_capacity/max_Capacity_N


stoi_n_data = [0];%:0.01:(1-N_stoichiometry_diff_0_100_SOC)];
stoi_p_data = [0.38];%:0.01:0.4];%(1-P_stoichiometry_diff_0_100_SOC)];
close all
% stoi_n_data = [0];%:0.01:(1-N_stoichiometry_diff_0_100_SOC)];
% stoi_p_data = [0:0.01:0.2];%:0.01:0.4];%(1-P_stoichiometry_diff_0_100_SOC)];
data=[];

for k=1:length(stoi_n_data)
    for m=1:length(stoi_p_data)

% 0% SOC
nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*stoi_n_data(k) + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*(stoi_p_data(m)+P_stoichiometry_diff_0_100_SOC); %match 0.000, 0.4835

% 100% SOC
nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*(stoi_n_data(k)+N_stoichiometry_diff_0_100_SOC) + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*stoi_p_data(m);%match 0.5647, 0.0115

clear ndata pdata OCV OCV_upp OCV_low

ndata = stoi_n_data(k):0.005:(stoi_n_data(k)+N_stoichiometry_diff_0_100_SOC);
nLis_avg = (nLis_0 + nLis_100)/2;%0.1796;
for i=1:length(ndata)
            
        theta_n = ndata(i);
        theta_p = (nLis_avg - p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*theta_n)/(p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max);
        pdata(i) = theta_p;
        
        theta_p
 OCPn(i) = refPotentialAnode(p,theta_n);
 OCPp(i) = refPotentialCathode(p,theta_p);
 
%  pocp = OCPp(i);
%  
%  pocp
 
 OCV(i) = OCPp(i)-OCPn(i);
 
%  figure(345)
%  plot3(ndata(i),pdata(i),OCV(i),'o')
%  hold on
%  xlabel('N')
%  ylabel('P')
   
end

%create SOC vector
SOC_vector = [0:1/(length(ndata)-1):1];


Vmean3 = interp1(SOC_data,Vmean2,SOC_vector);
% OCPn1 = interp1(SOC_vector,OCPn,SOC_data);
% OCPp1 = interp1(SOC_vector,OCPp,SOC_data);

norm_error = norm(abs(OCV-Vmean3),2);
figure(12)
plot(SOC_vector,OCV,'linewidth',2)
hold on
grid on
plot(SOC_vector,Vmean3,'r','linewidth',2)
xlabel('SOC')
ylabel('Open Circuit Voltage [V]')
legend('Est','Meas')

error = OCV - Vmean3;

data = [data; norm_error stoi_n_data(k) stoi_p_data(m) nLis_0 nLis_100];

    end
end
% plot(Vmean3)
% hold on
figure(36)
plot(SOC_vector,error,'b','linewidth',2)
grid on
xlabel('SOC')
ylabel('Error [V]')
title('Error = OCV-measured - OCV-fitted')


% figure(13)
% plot(SOC_data,OCV1 - error,'linewidth',2)
% hold on
% grid on
% plot(SOC_data,Vmean2,'r','linewidth',2)
% xlabel('SOC')
% ylabel('Open Circuit Voltage [V]')
% legend('Est','Meas')


figure(4)
subplot(2,1,1)
plot(pdata,OCPp,'linewidth',2)
grid on
xlabel('Normalized Concentration')
legend('OCP-P [V]')
subplot(2,1,2)
plot(ndata,OCPn,'linewidth',2)
xlabel('Normalized Concentration')
grid on
legend('OCP-N [V]')

% clear error_N error_P
% 
% error_N = error(1,995);
% error_P = zeros(1,995);
% % 
% % error_N(1:250) = error(1:250);
% % error_N(521:995) = error(521:995);
% % error_P(251:520) = error(1,251:520);;
% % 
% OCP_P_new = OCPp1+error_P;
% OCP_N_new = OCPn1+error_N;
% % 
% % 
% % 
% figure(5)
% subplot(2,1,1)
% plot(diff(OCP_P_new),'linewidth',2)
% grid on
% legend('diff-OCP-P [V]')
% subplot(2,1,2)
% plot(diff(OCP_N_new),'linewidth',2)
% grid on
% legend('diff-OCP-N [V]')
% ylim([-10e-4 5e-4])
% 
% 
% figure(6)
% subplot(2,1,1)
% plot(SOC_data,(OCP_P_new),'linewidth',2)
% grid on
% legend('diff-OCP-P [V]')
% subplot(2,1,2)
% plot(SOC_data,(OCP_N_new),'linewidth',2)
% grid on
% legend('diff-OCP-N [V]')


% 
% figure(3)
% plotyy(SOC_data,OCPp1,SOC_data,OCPn1)