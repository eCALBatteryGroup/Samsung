%%ACKLEY	Ackley's test function for numerical optimization
%
%	Known global minimum is the origin (in any dimension)
%
%	The function is vecotrized and supports any number of dimensions.
%
% Author(s):	Soren Ebbesen, 14-Sep-2011
%				sebbesen@idsc.mavt.ethz.ch
%
% Version:	1.1 (19-Jun-2013)
%
% Institute for Dynamic Systems and Control, Department of Mechanical and
% Process Engineering, ETH Zurich
%
% This Source Code Form is subject to the terms of the Mozilla Public License,
% v. 2.0. If a copy of the MPL was not distributed with this file, You can
% obtain one at http://mozilla.org/MPL/2.0/.


function f = DFN_error(x)

% % Dimensions
% n = size(x,2);
% 
% % Ackley's function (non-linear test function)
% f = 20 + exp(1) ...
%    -20*exp(-0.2*sqrt((1/n).*sum(x.^2,2))) ...
%    -exp((1/n).*sum(cos(2*pi*x),2));
disp('=================Parameters=================')
p.sig_n = x(1) %100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = x(2) %100;    % Conductivity of solid in pos. electrode, [1/Ohms*m]
% p.R_f_n = x(3); % Resistivity of SEI layer, [Ohms*m^2]
% p.R_s_n = x(4);%5.00E-7;%5.00E-06; %3.596e-6;   % Radius of solid particles in negative electrode [m]
% p.R_s_p = x(5);%5.00E-7; %1.637e-7;   % Radius of solid particles in positive electrode [m]
% p.R_f_p = x(6);%30.00E-03; %0;          % Resistivity of SEI layer, [Ohms*m^2]

run dfn_scott_testing_Satadru

load Data_Satadru_short_UDDS_100sec.mat

% size1 = size(out.volt);
% size2 = size(data(:,2));

% f = [out.time out.volt];
% figure(1)
% plot((out.volt - data(:,2)),'col',rand(1,3))
% hold on
% grid on
% 
% figure(2)
% plot(out.time,(out.volt),'b',out.time,data(:,2),'r','linewidth',2)
% grid on
% hold on
% plot((data(:,2)),'r')

f = norm((out.volt-data(:,2)),2);

%f = norm((Volt-volt_exp),2);


% figure
% plot(x(1),'*')
% figure
% plot(x(2),'*')
% disp('Design Variable Value');
% sig_n_value = x(1)
% sig_p_value = x(2)
%f = rms((out.volt - data(:,2)));