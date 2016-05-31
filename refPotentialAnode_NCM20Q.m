%% Reference Potential for Neg Electrode: Unref(theta_n)
%   Samsung INR1865020Q
%
%   Created May 21, 2016 by Scott Moura

function [Uref,varargout] = refPotentialAnode_NCM20Q(p,theta)

if(~isreal(theta))
    beep;
    error('dfn:err','Complex theta_n');
%     pause;
end

x = theta;

% Gaussian Mixture Model
%      General model Gauss6:
%      fittedmodel(x) = 
%               a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + 
%               a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + 
%               a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2)
%     Coefficients (with 95% confidence bounds):
a1 =   1.131e+07;  %(-1.582e+10, 1.584e+10)
b1 =     -0.1451;  %(-12.03, 11.74)
c1 =     0.03473;  %(-1.422, 1.491)
a2 =       1.032;  %(-7.247, 9.312)
b2 =    -0.02927;  %(-0.3278, 0.2693)
c2 =     0.04502;  %(-0.09087, 0.1809)
a3 =     0.06713;  %(-0.6805, 0.8147)
b3 =     0.06537;  %(-0.02014, 0.1509)
c3 =     0.03598;  %(-0.05277, 0.1247)
a4 =     0.08098;  %(-0.6746, 0.8366)
b4 =     0.07157;  %(-0.4527, 0.5959)
c4 =     0.07171;  %(-0.1449, 0.2883)
a5 =      0.1725;  %(0.1679, 0.1771)
b5 =      0.2342;  %(0.2178, 0.2505)
c5 =       0.275;  %(0.26, 0.2901)
a6 =   -0.002927;  %(-7.518e+05, 7.518e+05)
b6 =      0.1974;  %(-6822, 6822)
c6 =     0.00165;  %(-9.24e+04, 9.24e+04)

% Anode OCP
Uref = a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2) + ...
    a3*exp(-((x-b3)/c3).^2) + a4*exp(-((x-b4)/c4).^2) + ...
    a5*exp(-((x-b5)/c5).^2) + a6*exp(-((x-b6)/c6).^2);



% Gradient of OCP wrt c_s!!!
if(nargout >= 2)
%      df(x) = \Sum_i -ai*exp(-((x-bi)/ci)^2)*2*(x-bi)/ci^2

dUref_dtheta = -a1*exp(-((x-b1)/c1).^2)*2.*(x-b1)/c1^2 ...
    -a2*exp(-((x-b2)/c2).^2)*2.*(x-b2)/c2^2 ...
    -a3*exp(-((x-b3)/c3).^2)*2.*(x-b3)/c3^2 ...
    -a4*exp(-((x-b4)/c4).^2)*2.*(x-b4)/c4^2 ...
    -a5*exp(-((x-b5)/c5).^2)*2.*(x-b5)/c5^2 ...
    -a6*exp(-((x-b6)/c6).^2)*2.*(x-b6)/c6^2;
varargout{1} = dUref_dtheta / p.c_s_n_max;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end

