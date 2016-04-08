% Initial guess

Ar = 0.0102041;


Ar = 0.045;
esn = 0.7215;
esp = 0.6516;

Qn = esn*p.L_n*p.Faraday*0.048965*p.c_s_n_max*(0.8 - 0.001)/3600
Qp = esp*p.L_p*p.Faraday*0.044688*p.c_s_p_max*(0.95 - 0.359749)/3600


Qn_new = esn*p.L_n*p.Faraday*0.048965*p.c_s_n_max*1/3600
Qp_new = esp*p.L_p*p.Faraday*0.044688*p.c_s_p_max*1/3600




nLis_1 = esn*p.L_n*Ar*p.c_s_n_max*(0.790813 - 0.3);

clear ndata pdata OCV OCV_upp OCV_low

% ndata = 0:0.05:0.9;
% pdata = 0:0.05:1;
% 
% 
% for i=1:length(ndata)
%     for j=1:length(pdata)
%         
%         theta_n = ndata(i);
%         theta_p = pdata(j);
%         
%         
%  OCPn = refPotentialAnode(p,theta_n);
%  OCPp = refPotentialCathode(p,theta_p);
%  
%  OCV(i,j) = OCPp-OCPn;
%  
%     end
% end
% 
% mesh(ndata,pdata,OCV')
% xlabel('N')
% ylabel('P')
% zlim([2 4.2])
% ylim([0.3 1])


ndata = 0:0.05:1;
pdata = 0:0.05:0.55;
OCV_upp=[];
OCV_low=[];



for i=1:length(ndata)
    for j=1:length(pdata)
        
        theta_n = ndata(i);
        theta_p = pdata(j);
        
        
 OCPn = refPotentialAnode(p,theta_n);
 OCPp = refPotentialCathode(p,theta_p);
 
 OCV(i,j) = OCPp-OCPn;
 
 if OCV(i,j)<=4.22 && OCV(i,j)>=4.18
     
     OCV_upp = [OCV_upp; OCV(i,j) theta_n theta_p];
     
 end
 
 if OCV(i,j)<=2.3 && OCV(i,j)>=2
     
     OCV_low = [OCV_low; OCV(i,j) theta_n theta_p];
     
 end
 
    end
    
    plot(pdata, OCV(i,:))
    hold on
end

plot([0 1],[4.2 4.2], '--k')
plot([0 1],[2.0 2.0], '--k')

% 0% SOC
nLis_0 = esn*p.L_n*Ar*p.c_s_n_max*0.001 + esp*p.L_p*Ar*p.c_s_p_max*0.95

% 100% SOC
nLis_100 = esn*p.L_n*Ar*p.c_s_n_max*0.8 + esp*p.L_p*Ar*p.c_s_p_max*0.35

clear ndata pdata OCV OCV_upp OCV_low

ndata = 0:0.005:0.8;
nLis_avg = 0.1672;%0.1796;
for i=1:length(ndata)
            
        theta_n = ndata(i);
        theta_p = (nLis_avg - esn*p.L_n*Ar*p.c_s_n_max*theta_n)/(esp*p.L_p*Ar*p.c_s_p_max);
        pdata(i) = theta_p;
        
 OCPn = refPotentialAnode(p,theta_n);
 OCPp = refPotentialCathode(p,theta_p);
 
 OCV(i) = OCPp-OCPn;
   
end
figure(12)
plot(ndata,OCV,'linewidth',2)
xlabel('Norm. Li concentration in Anode')
ylabel('Open Circuit Voltage [V]')



%% From Tanim et al 2015 data

Qn = 0.662*40e-4*96487*1020.41*31.08e-3*(0.8 - 0.001)/3600
Qp = 0.58*36.55e-4*96487*1020.41*51.83e-3*(0.95 - 0.359749)/3600


Qn_new = 0.662*40e-4*96487*1020.41*31.08e-3*1/3600
Qp_new = 0.58*36.55e-4*96487*1020.41*51.83e-3*1/3600


