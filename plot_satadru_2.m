sum=0;
SOC_CC(1)=0.7674;
SOC_CC2(1)=0.7674;
for k = 1:(NT-1)
   SOC_cathode(k+1) = mean(c_avg_p(:,k+1)) / p.c_s_p_max; 
   
   SOC_CC(k+1)=SOC_CC(k)-I(k)*p.Area*p.delta_t/(OneC*3600);
   SOC_CC2(k+1)=SOC_CC2(k)-I(k)*p.Area*p.delta_t/(2.89*3600);
%    SOC_CC(k+1)=SOC_CC(k)-I(k)*0.048965*p.delta_t/(3.62*3600);
   
%    sum = sum + I(k);
%     
%    jn_mean(k+1)=mean(jn(:,k+1));
%    
%    i_en_mean(k+1)=mean(i_en(:,k+1));
   
   SOC_anode(k+1) = mean(c_avg_n(:,k+1))/((0.575)*p.c_s_n_max);
end


% figure(34)
% clf
% plot(t,SOC_CC,'linewidth',2)
% legend('SOC - Coulomb Counting')
% 
% 
% figure(36)
% clf
% plot(t,i_en_mean,'linewidth',2)
% hold on
% plot(t,I,'r','linewidth',2)
% legend('I from jn','actual I')
% 
% figure(35)
% clf
% plot(t,jn_mean*(p.epsilon_s_n*p.L_n*p.Area*p.Faraday),'linewidth',2)
% hold on
% plot(t,I,'r','linewidth',2)
% legend('I from jn','actual I')


% figure(1)
% clf
% subplot(211)
% plot(t,(I*p.Area)/OneC,'linewidth',2)
% legend('I [C-rate]')
% title('Current')
% xlim([0 t(end)])
% subplot(212)
% plot(t,Volt,'linewidth',2)
% legend('V [V]')
% ylim([3.5 3.9])
% xlim([0 t(end)])
% title('Voltage')
% xlabel('Time [s]')


% figure(2)
% plotyy(t(2:end),SOC(2:end),t(2:end),SOC_cathode(2:end))
% hold on
% %plot(t,SOC_CC,'--r','linewidth',2)
% 
% ylabel('Cell SOC [0-1]')
% title('SOC [0-1]')
% xlabel('Time [s]')

figure(21)
plot(t(2:end),SOC_anode(2:end),'linewidth',2)
hold on
plot(t,SOC_CC,'--k','linewidth',2)
plot(t,SOC_CC2,'--r','linewidth',2)
ylabel('Cell SOC [0-1]')
title('SOC [0-1]')
xlabel('Time [s]')
legend('DFN: With Equivalent Capacity 2.89 Ah','CC: With Rated Capacity 2.08 Ah','CC: With Capacity 2.89 Ah')
%ylim([0.453 0.5518])
figure(22)
plot(t,SOC_CC-SOC_anode,'linewidth',2)



% figure(3)
% subplot(4,1,1)
% plot(t,c_ss_n,'linewidth',2)
% hold on
% ylabel('Li con. [mol/m^3]')
% title('Negative Electrode: Surface Concentration @ each node')
% subplot(4,1,2)
% plot(t,c_ss_p,'linewidth',2)
% title('Positive Electrode: Surface Concentration @ each node')
% ylabel('Li con. [mol/m^3]')
% %xlabel('Time [s]')
% subplot(413)
% plot(t,c_e_0p/1e3,'linewidth',2)
% hold on
% %plot(t,0.15*ones(size(t)),'k--')
% legend('ce0p [kmol/m^3]')
% xlim([0 t(end)])
% title('Li concentration in electrolyte @ boundary')
% subplot(414)
% plot(t,c_e(:,:),'linewidth',2)
% hold on
% %plot(t,0*ones(size(t)),'k--')
% %legend('eta')
% % ylim([-0.15 0.2])
% xlim([0 t(end)])
% ylabel('Li con. [mol/m^3]')
% xlabel('Time [s]')
% title('Li concentration in electrolyte @ all nodes')