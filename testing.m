load Vmean2_Co20_end2p6.mat
close all
pointdata=[0.02:0.005:0.05];
pointdata=[0.06:0.005:0.2];
for k = 1:length(pointdata)

point = 0.02;%pointdata(k);
% 0% SOC
nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*point + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.4835 %match 0.000, 0.4835

% 100% SOC
nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*(0.5647+point) + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.0115 %match 0.5647, 0.0115

clear ndata pdata OCV OCV_upp OCV_low OCPn

ndata = point:0.005:0.5647+point;
%ndata = 0.0:0.005:1;

nLis_avg = pointdata(k);%(nLis_0 + nLis_100)/2;%0.1796;
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
% figure(12)
% plot(SOC_vector,OCPn,'r','linewidth',2)
% hold on
% grid on
% plot(SOC_data,OCP_Gr,'b','linewidth',2)

figure(13)
plot(SOC_data,OCP_Gr,'col',rand(3,1),'linewidth',2)
hold on
grid on
%plot(SOC_data,Vmean2,'r','linewidth',2)
xlabel('Norm. Li concentration in Anode')
ylabel('Open Circuit Voltage [V]')
legend('Est','Meas')



figure(14)
plot(SOC_data,OCP_Gr + Vmean2,'col',rand(3,1),'linewidth',2)
hold on
grid on
% plot(SOC_data,Vmean2,'r','linewidth',2)
% plot(SOC_data,-OCP_Gr,'k','linewidth',2)
xlabel('Norm. Li concentration in Anode')
ylabel('OCP-P')


% figure(15)
% plot(SOC_data,OCP_Gr + Vmean2,'col',rand(3,1),'linewidth',2)
% hold on
% plot(SOC_data,OCP_Gr,'col',rand(3,1),'linewidth',2)
% plot(SOC_data,Vmean2,'r','linewidth',2)
% 
% plotyy(SOC_data,OCP_Gr + Vmean2,SOC_data,-OCP_Gr)
% hold on
% plot(SOC_data,Vmean2,'r','linewidth',2)




end