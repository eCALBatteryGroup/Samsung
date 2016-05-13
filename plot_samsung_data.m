%% Plot Experimental Data collected for Samsung
%   Created May 12, 2016 by Scott Moura

close all;
fs = 16;

figure(1); clf;



load('data/Int_Obs/1C_data_Oct_26_2015_05_sample')

subplot(3,2,1)
plot(time_exp,current_exp,'LineWidth',2)
xlim([time_exp(1), time_exp(end)])
ylim([-4,0])
set(gca,'FontSize',fs)

title('\bfCurrent [A]','FontSize',fs)

subplot(3,2,2)
plot(time_exp,volt_exp,'LineWidth',2)
xlim([time_exp(1), time_exp(end)])
set(gca,'FontSize',fs)

title('\bfVoltage [V]','FontSize',fs)

clear time_exp current_exp volt_exp



load('data/Int_Obs/dfn_5c')

subplot(3,2,3)
plot(time_exp,current_exp,'LineWidth',2)
xlim([time_exp(1), time_exp(end)])
set(gca,'FontSize',fs)

subplot(3,2,4)
plot(time_exp,volt_exp,'LineWidth',2)
xlim([time_exp(1), time_exp(end)])
set(gca,'FontSize',fs)

clear time_exp current_exp volt_exp



load('data/Int_Obs/UDDS_data_Oct_26_2015_Sample_05sec');

subplot(3,2,5)
plot(time_exp,current_exp,'LineWidth',2)
xlim([time_exp(1), time_exp(end)])
set(gca,'FontSize',fs)
xlabel('Time [sec]','FontSize',fs)

subplot(3,2,6)
plot(time_exp,volt_exp,'LineWidth',2)
xlim([time_exp(1), time_exp(end)])
set(gca,'FontSize',fs)
xlabel('Time [sec]','FontSize',fs)