clc
close all

x_data = [0.001:0.01:1];

for i=1:length(x_data)
    
OCP_NMC(i)= -10.72*x_data(i)^4 + 23.88*x_data(i)^3 -16.77*x_data(i)^2 + 2.595*x_data(i) +4.563;    
    %OCPp(i) = interp1((1-Norm_concentration),OCP_P,x_data(i));
 
end

figure(1)
plot(x_data,OCP_NMC,'linewidth',2)
grid on
hold on
%plot(x_data,OCPp,'r','linewidth',2)
legend('Literature','Samsung')
title('OCP: NMC Electrode')
xlabel('Normized Li Concentration [0-1]')
ylabel('Potential [V]')
%ylim([3.5 4.5])



x_data2 = [0.001:0.01:1];

for i=1:length(x_data2)
    
OCP_G(i)= 0.1493 + 0.8493*exp(-61.79*x_data2(i))+0.3824*exp(-665.8*x_data2(i))- exp(39.42*x_data2(i)-41.92)-0.03131*atan(25.59*x_data2(i)-4.099)-0.009434*atan(32.49*x_data2(i)-15.74);  
  %OCPn(i) = interp1(Norm_concentration,OCP_N,x_data2(i)); 


end

figure(2)
plot(x_data2,OCP_G,'linewidth',2)
grid on
hold on
%plot(x_data2,OCPn,'r','linewidth',2)
legend('Literature','Samsung')
title('OCP: Graphite Electrode')
xlabel('Normized Li Concentration [0-1]')
ylabel('Potential [V]')
% ylim([3.5 4.5])