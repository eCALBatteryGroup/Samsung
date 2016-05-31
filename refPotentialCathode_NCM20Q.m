%% Reference Potential for Positive Electrode: Upref(theta_p)
%   Samsung INR1865020Q
%
%   Created May 21, 2016 by Scott Moura

function [Uref,varargout] = refPotentialCathode_NCM20Q(p,theta)

if(~isreal(theta))
    beep;
    error('dfn:err','Complex theta_n');
%     pause;
end

x = theta;

% Polynomail Model
%      Linear model Poly4:
%      fittedmodel(x) = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5
%      Coefficients (with 95% confidence bounds):
p1 = -10.715239318554429;  %(-10.72, -10.71)
p2 = 23.869645625119176;   %(23.87, 23.87)
p3 = -16.762450136934714;  %(-16.76, -16.76)
p4 = 2.593038841632656;    %(2.592, 2.594)
p5 = 4.563097845038069;    %(4.563, 4.563)

% Anode OCP
Uref = p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x + p5;



% Gradient of OCP wrt c_s!!!
if(nargout >= 2)
%      df(x) = 4*p1*x.^3 + 3*p2*x.^2 + 2*p3*x + p4

dUref_dtheta = 4*p1*x.^3 + 3*p2*x.^2 + 2*p3*x + p4;
varargout{1} = dUref_dtheta / p.c_s_p_max;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end

