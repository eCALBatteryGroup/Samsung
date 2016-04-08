load Vmean2_Co20_end2p6.mat
close all
% figure(3)
% plot(SOC_data,Vmean2)

pp = 0.3;

point = 0.05;
% 0% SOC
nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*point + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*(pp+P_stoichiometry_diff_0_100_SOC) %match 0.000, 0.4835

% 100% SOC
nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*(N_stoichiometry_diff_0_100_SOC+point) + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*(pp) %match 0.5647, 0.0115

clear ndata pdata OCV OCV_upp OCV_low OCPn

ndata = point:0.005:0.5647+point;
%ndata = 0.0:0.005:1;

nLis_avg = (nLis_0 + nLis_100)/2;%0.1796;
for i=1:length(ndata)
            
        theta_n = ndata(i);
        theta_p = (nLis_avg - p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*theta_n)/(p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max);
        pdata(i) = theta_p;
        
 OCPn(i) = refPotentialAnode(p,theta_n);
%  OCPp = refPotentialCathode(p,theta_p);
%  
%  OCV(i) = OCPp-OCPn;
   
end

% plot(ndata,OCPn,'r')
% hold on
clear SOC_vector
%create SOC vector
SOC_vector = [0.001:1/length(OCPn):0.999];

OCP_Gr = interp1(SOC_vector,OCPn,SOC_data);
figure(12)
plot(SOC_vector,OCPn,'r','linewidth',2)
hold on
grid on
plot(SOC_data,OCP_Gr,'b','linewidth',2)

% figure(13)
% plot(SOC_data,OCP_Gr,'b','linewidth',2)
% hold on
% grid on
% plot(SOC_data,Vmean2,'r','linewidth',2)
% xlabel('Norm. Li concentration in Anode')
% ylabel('Open Circuit Voltage [V]')
% legend('Est','Meas')



figure(14)
plot(SOC_data,OCP_Gr + Vmean2,'linewidth',2)
hold on
grid on
% plot(SOC_data,Vmean2,'r','linewidth',2)
% plot(SOC_data,-OCP_Gr,'k','linewidth',2)
xlabel('Norm. Li concentration in Anode')
ylabel('Open Circuit Voltage [V]')


figure(15)
plotyy(SOC_data,OCP_Gr + Vmean2,SOC_data,OCP_Gr)

figure(16)
plot(diff(OCP_Gr + Vmean2),'linewidth',2)
hold on
grid on
ylim([-0.5e-3 2e-3])