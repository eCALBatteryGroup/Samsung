%% Varphi Generator
%   Script to generate c_ss_n = varphi(V,I) lookup table
%   Created Jun 26, 2015 by Scott Moura

clear;

%% Model Construction
% Electrochemical Model Parameters
% run params_bosch
run params_FePO4_ACC15

V_vec = linspace(2.0,3.6,501);
I_vec = linspace(-20,20,51);
% I_vec = 0;

%% Numerically determine c_ss_n from (V,I) pair

RTaF=(p.R*p.T_amp)/(p.alph*p.Faraday);

varphi = nan*zeros(length(V_vec),length(I_vec));

for idx1 = 1:length(V_vec);
    
    for idx2 = 1:length(I_vec);
        
        V = V_vec(idx1);
        I = I_vec(idx2);
        
        % Algorithm params
        maxiters = 50;
        x = zeros(maxiters,1);
        f = nan*ones(maxiters,1);
        tol = 1e-5;

        % Initial Guesses
        x_low = 0.0 * p.c_s_n_max;
        x_high = 1.0 * p.c_s_n_max;
        x(1) = 0.5 * p.c_s_n_max;
        
        % Iterate Bisection Algorithm
        for idx = 1:maxiters

            % Stochiometric SOC
            theta_n = x(idx)/p.c_s_n_max;
            theta_p = (p.n_Li_s-p.epsilon_s_n*p.L_n*p.Area*x(idx))/(p.c_s_p_max*p.epsilon_s_p*p.L_p*p.Area);

            % OCP Functions
            OCPn = refPotentialAnode(p,theta_n);
            OCPp = refPotentialCathode(p,theta_p);

            % Exchange Current Density
            c_e0 = p.c_e * ones(p.Nx,1);     % Fixed electrolyte concentration [mol/m^3]
            [i_0n,i_0p] = exch_cur_dens(p,x(idx),theta_p*p.c_s_p_max,c_e0);

            % Compute output function
            f(idx) = RTaF * asinh(-I / (2*p.a_s_p*p.Area*p.L_p*i_0p(end))) ...
                    -RTaF * asinh(I / (2*p.a_s_n*p.Area*p.L_n*i_0n(1))) ...
                    + OCPp - OCPn ...
                    - (p.R_f_n/(p.a_s_n*p.L_n*p.Area) + p.R_f_p/(p.a_s_p*p.L_p*p.Area))*I - V;

            % Bisection
            if(abs(f(idx)) <= tol)
                break;
            elseif(f(idx) >= 0)
                x_high = x(idx);
            else
                x_low = x(idx);
            end

            % Bisection
            x(idx+1) = (x_high + x_low)/2;
            x(idx+1)/p.c_s_n_max;

        end

        % Output conveged csp0
%         theta_n = theta_n
%         theta_p = (p.n_Li_s-p.epsilon_s_n*p.L_n*x(idx))/(p.c_s_p_max*p.epsilon_s_p*p.L_p)
        varphi(idx1,idx2) = x(idx);
        
    end

end

%%
fs = 16;

figure(1); clf;
surf(V_vec,I_vec,varphi'/1e3);
xlim([2.0 3.6]);
ylim([-20 20]);
xlabel('Voltage [V]','FontSize',fs);
ylabel('Current [A]','FontSize',fs);
% zlabel('$$\varphi(V,I) [kmol/m^3]$$','Interpreter','latex','FontSize',fs);
set(gca,'FontSize',fs);

