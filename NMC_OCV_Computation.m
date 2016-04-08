clc
clear all
close all
params_NMC_Samsung

%% Step 1: Covert SAMASUNG data to uniform "concentration vs OCP" functions
% load OCP_data_SAMSUNG.mat
% 
% Norm_concentration = [0:0.01:1];
% OCP_P = interp1(SOC_Cathode,OCP_Cathode,Norm_concentration);
% OCP_N = interp1(SOC_Anode,OCP_Anode,Norm_concentration);
% 
% figure(1)
% plot(SOC_Cathode,OCP_Cathode)
% hold on
% grid on
% plot(Norm_concentration,OCP_P,'r')
% 
% figure(2)
% plot(SOC_Anode,OCP_Anode)
% hold on
% grid on
% plot(Norm_concentration,OCP_N,'r')

%% Step 2: Maximum concentration calculation from SAMSUNG OCP data

load OCP_data_NMC_SAMSUNG.mat

cell_capacity = 2.11;%Ah
max_Capacity_P = 4.34;%Ah
max_Capacity_N = 3.63;%Ah

max_concentration_P = (max_Capacity_P*3600)/(p.epsilon_s_p*p.L_p*p.Area*p.Faraday);%mol/m3
max_concentration_N = (max_Capacity_N*3600)/(p.epsilon_s_n*p.L_n*p.Area*p.Faraday);%mol/m3

P_stoichiometry_diff_0_100_SOC = cell_capacity/max_Capacity_P
N_stoichiometry_diff_0_100_SOC = cell_capacity/max_Capacity_N


%% Step 3:Measured OCV Curve

load Vmean2_Co20_end2p6.mat
Vmean2 = Vmean2(5:end);
SOC_data = [0:1/(length(Vmean2)-1):1];

figure(3)
plot(SOC_data,Vmean2)

%% Step 4: Find stoichiometry points

 
clear ndata pdata OCV OCV_upp OCV_low
ndata = 0:0.01:1;
pdata = 0:0.01:1;
OCV_upp=[];
OCV_low=[];



for i=1:length(ndata)
    for j=1:length(pdata)
        
        theta_n = ndata(i);
        theta_p = pdata(j);
        
        
 OCPn = refPotentialAnode(p,theta_n);
 OCPp = refPotentialCathode(p,theta_p);
 
%  OCPp = interp1(Norm_concentration,OCP_P,theta_p);
%  OCPn = interp1(Norm_concentration,OCP_N,theta_n);
 
 OCV(i,j) = OCPp-OCPn;
 
 if OCV(i,j)<=4.22 && OCV(i,j)>=4.18
     
     OCV_upp = [OCV_upp; OCV(i,j) theta_n theta_p];
     
 end
 
 if OCV(i,j)<=2.8 && OCV(i,j)>=2
     
     OCV_low = [OCV_low; OCV(i,j) theta_n theta_p];
     
 end
 
    end
    
%     plot(pdata, OCV(i,:))
%     hold on
end

mesh(ndata,pdata,OCV')
xlabel('Norm. c-ss-N [0-1]')
ylabel('Norm. c-ss-P [0-1]')
zlabel('OCV = OCP_P - OCP_N')


%% Using stoichiometry points, find n-Li-s and plot OCV
% 0% SOC
nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*0.000 + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.4835 %match 0.000, 0.4835

% 100% SOC
nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*0.5647 + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.0115%match 0.5647, 0.0115

clear ndata pdata OCV OCV_upp OCV_low

ndata = 0:0.005:0.59;
nLis_avg = (nLis_0 + nLis_100)/2;%0.1796;
for i=1:length(ndata)
            
        theta_n = ndata(i);
        theta_p = (nLis_avg - p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*theta_n)/(p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max);
        pdata(i) = theta_p;
        
 OCPn = refPotentialAnode(p,theta_n);
 OCPp = refPotentialCathode(p,theta_p);
 
 OCV(i) = OCPp-OCPn;
   
end

%create SOC vector
SOC_vector = [0:1/(length(ndata)-1):1];

figure(12)
plot(SOC_vector,OCV,'linewidth',2)
hold on
grid on
plot(SOC_data,Vmean2,'r','linewidth',2)
xlabel('Norm. Li concentration in Anode')
ylabel('Open Circuit Voltage [V]')
legend('Est','Meas')


%% To find a n-Li-s that gives a good OCV fit

load Vmean2_Co20_end2p6.mat


%SOC_data = [0:1/(length(Vmean2)-1):1];


clear ndata SOC_vector
figure(12)
plot(SOC_data,Vmean2,'r','linewidth',2)
hold on
grid on
xlabel('Norm. Li concentration in Anode')
ylabel('Open Circuit Voltage [V]')
%legend('Est','Meas')
ndata = [0:0.005:0.59];
nLidata = [0.1111];
for j=1:length(nLidata)
nLis_avg = nLidata(j);%0.1796;
for i=1:length(ndata)
            
        theta_n = ndata(i);
        theta_p = (nLis_avg - p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*theta_n)/(p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max);
        pdata(i) = theta_p;
        
 OCPn = refPotentialAnode(p,theta_n);
 OCPp = refPotentialCathode(p,theta_p);
 
 OCV(i) = OCPp-OCPn;
   
end

%create SOC vector
SOC_vector = [0:1/(length(ndata)-1):1];

figure(12)
plot(SOC_vector,OCV,'col',rand(3,1),'linewidth',2)
hold on
grid on

end

%% Discard Samsung OCP data, use literature OCP-n and subtract from OCV to get OCP-p


load Vmean2_Co20_end2p6.mat

figure(3)
plot(SOC_data,Vmean2)

point = 0.07;
% 0% SOC
nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*point + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.4835 %match 0.000, 0.4835

% 100% SOC
nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*(0.5647+point) + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.0115 %match 0.5647, 0.0115

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

plot(ndata,OCPn,'r')
hold on
clear SOC_vector
%create SOC vector
SOC_vector = [0.001:1/length(OCPn):0.999];

OCP_Gr = interp1(SOC_vector,OCPn,SOC_data);
figure(12)
plot(SOC_vector,OCPn,'r','linewidth',2)
hold on
grid on
plot(SOC_data,OCP_Gr,'b','linewidth',2)

figure(13)
plot(SOC_data,OCP_Gr,'b','linewidth',2)
hold on
grid on
plot(SOC_data,Vmean2,'r','linewidth',2)
xlabel('Norm. Li concentration in Anode')
ylabel('Open Circuit Voltage [V]')
legend('Est','Meas')



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
%% Discard Samsung OCP data, use literature OCP-n and subtract from OCV to get OCP-p
% 
% 
% load Vmean2_Co20.mat
% Vmean2(1)=2;
% SOC_data = [0.001:0.001:0.999];
% 
% figure(3)
% plot(SOC_data,Vmean2)
% 
% % 0% SOC
% nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*0.000 + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.4835 %match 0.000, 0.4835
% 
% % 100% SOC
% nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*0.5647 + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.0115 %match 0.5647, 0.0115
% 
% clear ndata pdata OCV OCV_upp OCV_low
% 
% ndata = 0:0.005:0.5647;
% pdata = 0.0115:0.005:0.4835;
% 
% nLis_avg = (nLis_0 + nLis_100)/2;%0.1796;
% for i=1:length(pdata)
%             
% %         theta_n = ndata(i);
% %         theta_p = (nLis_avg - p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*theta_n)/(p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max);
% %         pdata(i) = theta_p;
%         
%  %OCPn(i) = refPotentialAnode(p,theta_n);
%  OCPp(i) = refPotentialCathode(p,pdata(i));
% %  
% %  OCV(i) = OCPp-OCPn;
%    
% end
% clear SOC_vector
% %create SOC vector
% SOC_vector = 1-[0.001:1/length(OCPp):0.999];
% 
% OCP_NMC = interp1(SOC_vector,OCPp,SOC_data);
% figure(12)
% plot(SOC_vector,OCPp,'r','linewidth',2)
% hold on
% grid on
% plot(SOC_data,OCP_NMC,'b','linewidth',2)
% 
% figure(13)
% plot(SOC_data,OCP_NMC,'b','linewidth',2)
% hold on
% grid on
% plot(SOC_data,Vmean2,'r','linewidth',2)
% xlabel('Norm. Li concentration in Anode')
% ylabel('Open Circuit Voltage [V]')
% legend('Est','Meas')
% 
% 
% 
% figure(14)
% plot(SOC_data,OCP_NMC - Vmean2,'linewidth',2)
% hold on
% grid on
% plot(SOC_data,Vmean2,'r','linewidth',2)
% %plot(SOC_data,-OCP_Gr,'k','linewidth',2)
% xlabel('Norm. Li concentration in Anode')
% ylabel('Open Circuit Voltage [V]')
















