p1data = [0:0.05:(1-0.5647)];
n1data = [0:0.05:(1-0.4724)];
for m=1:length(p1data)
    for k=1:length(n1data)
% 0% SOC
nLis_0 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*0.000 + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.4835 %match 0.000, 0.4835

% 100% SOC
nLis_100 = p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*0.5647 + p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max*0.0115%match 0.5647, 0.0115

clear ndata pdata OCV OCV_upp OCV_low

ndata = n1:0.005:0.5647+n1;
pdata = p1:0.005:0.4724+p1;
nLis_avg = (nLis_0 + nLis_100)/2;%0.1796;
for i=1:length(ndata)
    for j=1:length(pdata)
            
        theta_n = ndata(i);
        theta_p = pdata(i);
%         theta_p = (nLis_avg - p.epsilon_s_n*p.L_n*p.Area*p.c_s_n_max*theta_n)/(p.epsilon_s_p*p.L_p*p.Area*p.c_s_p_max);
%         pdata(i) = theta_p;
        
 OCPn(i) = refPotentialAnode(p,theta_n);
 OCPp(i) = refPotentialCathode(p,theta_p);
 
 OCV(i) = OCPp(i)-OCPn(i);
    end
end
clear SOC_vector
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

end
    end
