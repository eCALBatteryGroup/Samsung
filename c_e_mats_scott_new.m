%% Matrices for Li Diffusion in Electrolyte Phase, c_e(x,t)
%   Created May 21, 2016 by Scott Moura
%
%   A complete reboot of the original c_e_mats from July 2011
%   Updated to incorporate eletrode/separator BCs into differential eqns   

function [M1,M2,M3,M4,M5, C_ce] = c_e_mats_scott_new(p)

%% Lumped Coefficients
Del_xn = p.L_n * p.delta_x_n;
Del_xs = p.L_s * p.delta_x_s;
Del_xp = p.L_p * p.delta_x_p;

Del_xns = (Del_xn + Del_xs)/2;
Del_xsp = (Del_xs + Del_xp)/2;

%%% Matrices in nonlinear dynamics
%% M1*cx + M2*cz

%%% M1
% eqns for anode, n
M1n_l = diag(ones(p.Nxn-2,1),+1) - diag(ones(p.Nxn-2,1),-1);
M1n_c = zeros(p.Nxn-1,1);
M1n_c(end) = 1;

M1n = [M1n_l, M1n_c, zeros(p.Nxn-1, (p.Nxs-1 + 1 + p.Nxp-1))];
M1n = M1n*1/(2*Del_xn);


% eqns for ns interface
M1ns_c = [-1/(2*Del_xn), 1/(2*Del_xn)-1/(2*Del_xs), 1/(2*Del_xs)];
M1ns = [zeros(1,p.Nxn-2), M1ns_c, zeros(1,p.Nxs-2 + 1 + p.Nxp-1)];


% eqns for separator, s
M1s_c = diag(ones(p.Nxs-2,1),+1) - diag(ones(p.Nxs-2,1),-1);

M1s_l = zeros(p.Nxs-1,1);
M1s_l(1) = -1;
M1s_r = zeros(p.Nxs-1,1);
M1s_r(end) = 1;

M1s = [zeros(p.Nxs-1,p.Nxn-1), M1s_l, M1s_c, M1s_r, zeros(p.Nxs-1,p.Nxp-1)];
M1s = M1s*1/(2*Del_xs);


% eqns for sp interface
M1sp_c = [-1/(2*Del_xs), 1/(2*Del_xs)-1/(2*Del_xp), 1/(2*Del_xp)];
M1sp = [zeros(1,p.Nxn-1 + 1 + p.Nxs-2), M1sp_c, zeros(1,p.Nxp-2)];


% eqns for cathode, p
M1p_r = diag(ones(p.Nxp-2,1),+1) - diag(ones(p.Nxp-2,1),-1);
M1p_c = zeros(p.Nxp-1,1);
M1p_c(1) = -1;

M1p = [zeros(p.Nxp-1, (p.Nxn-1 + 1 + p.Nxs-1)), M1p_c, M1p_r];
M1p = M1p*1/(2*Del_xp);


% assemble submatrices
M1 = [M1n; M1ns; M1s; M1sp; M1p];
M1 = sparse(M1);


%%% M2
M2 = zeros(p.Nx-1,2);
M2(1,1) = -1/(2*Del_xn);
M2(end,end) = 1/(2*Del_xp);
M2 = sparse(M2);


%% M3*cx + M4*cz
%%% M3
% eqns for anode, n
M3n_l = -2*diag(ones(p.Nxn-1,1),0) + diag(ones(p.Nxn-2,1),+1) + diag(ones(p.Nxn-2,1),-1);
M3n_c = zeros(p.Nxn-1,1);
M3n_c(end) = 1;

M3n = [M3n_l, M3n_c, zeros(p.Nxn-1, (p.Nxs-1 + 1 + p.Nxp-1))];
M3n = M3n/(Del_xn^2);


% eqns for ns interface
M3ns_c = [1/(Del_xn*Del_xns), -(1/(Del_xn*Del_xns)+1/(Del_xs*Del_xns)), 1/(Del_xs*Del_xns)];
M3ns = [zeros(1,p.Nxn-2), M3ns_c, zeros(1,p.Nxs-2 + 1 + p.Nxp-1)];


% eqns for separator, s
M3s_c = -2*diag(ones(p.Nxs-1,1),0) + diag(ones(p.Nxs-2,1),+1) + diag(ones(p.Nxs-2,1),-1);

M3s_l = zeros(p.Nxs-1,1);
M3s_l(1) = 1;
M3s_r = zeros(p.Nxs-1,1);
M3s_r(end) = 1;

M3s = [zeros(p.Nxs-1,p.Nxn-1), M3s_l, M3s_c, M3s_r, zeros(p.Nxs-1,p.Nxp-1)];
M3s = M3s/(Del_xs^2);


% eqns for sp interface
M3sp_c = [1/(Del_xs*Del_xsp), -(1/(Del_xs*Del_xsp)+1/(Del_xp*Del_xsp)), 1/(Del_xp*Del_xsp)];
M3sp = [zeros(1,p.Nxn-1 + 1 + p.Nxs-2), M3sp_c, zeros(1,p.Nxp-2)];


% eqns for cathode, p
M3p_r = -2*diag(ones(p.Nxp-1,1),0) + diag(ones(p.Nxp-2,1),+1) + diag(ones(p.Nxp-2,1),-1);
M3p_c = zeros(p.Nxp-1,1);
M3p_c(1) = 1;

M3p = [zeros(p.Nxp-1, (p.Nxn-1 + 1 + p.Nxs-1)), M3p_c, M3p_r];
M3p = M3p/(Del_xp^2);


% assemble submatrices
M3 = [M3n; M3ns; M3s; M3sp; M3p];
M3 = sparse(M3);


%%% M4
M4 = zeros(p.Nx-1,2);
M4(1,1) = 1/(Del_xn^2);
M4(end,end) = 1/(Del_xp^2);
M4 = sparse(M4);


%% M5*j
M5n = (1-p.t_plus)*p.a_s_n/p.epsilon_e_n * speye(p.Nxn-1);
M5p = (1-p.t_plus)*p.a_s_p/p.epsilon_e_p * speye(p.Nxp-1);

rM5 = [p.Nxn-1; p.Nxs+1; p.Nxp-1];
cM5 = rM5';
M5 = blkdiagFast(rM5, cM5, M5n, zeros(p.Nxs+1), M5p);
M5 = sparse(M5);

%% Boundary Conditions
N1 = zeros(2,p.Nx-1);
N2 = zeros(2);

% BC1
N1(1,1) = +4;
N1(1,2) = -1;

N2(1,1) = -3;


% BC2
N1(2,end-1) = +1;
N1(2,end) = -4;

N2(2,2) = +3;

%%% SPARSE OUTPUT 
C_ce = sparse(-N2\N1);


