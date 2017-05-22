function solve_FD_spheres_variable_diffusivity(varargin)
%% Solves the PDE for Fickian diffusion with variable diffusivity within a radially symmetric sphere 
% The PDE is solved numerically using method of lines, the spherical finite
% difference discretization scheme defined and analyzed in the reference, 
% and any built-in MATLAB ODE solver. 
% Examples provided for variable diffusivity cases I-IV presented in the 
% reference listed below. The use of ode45 solver is shown in this code.  
% Any ODE solver could be selected. The reference used RADAU5 in Fortran.

%% --- Reference --- 
%  A. N. Ford Versypt & R. D. Braatz, Analysis of finite difference 
%  discretization schemes for diffusion in spheres with variable 
%  diffusivity, Computers & Chemical Engineering, 71 (2014) 241-252, 
%  https://doi.org/10.1016/j.compchemeng.2014.05.022

clf

%% --- Input ---
% Each case varies the functional form for the diffusivity alpha
if nargin==1
    % read DCASE from input to the function call in the command window
    DCASE = varargin{1};
else
    % specifiy DCASE
    % integer value 1, 2, 3, or 4 for Cases I, II, III, or IV described below and in the reference
    DCASE = 1; 
end
%% Diffusivity cases (defined computationally in FD_spheres_variable_diffusivity.m)
% DCASE = 1; 
% Case I: constant diffusivity with alpha = alpha0

% DCASE = 2;
% Case II: time-dependent diffusivity
%   alpha = alpha0 for t =< tau
%   and
%   alpha = alpha0*exp(k*(t - tau)) for t > tau

% DCASE = 3;
% Case III: spatially dependent diffusivity
%   alpha = alpha0*r^2

% DCASE = 4;
% CASE IV: concentration-dependent diffusivity
%   alpha = alpha0*exp(k-k*c)

%% Parameters
M = 100; %number of spatial intervals in Method of Lines
NR = M+1; % number of spatial discretization points along r
initial_condition = 1; % c(r,0) for 0 <= r < 1
boundary_condition = 0; % c(1,t) for t >= 0
c0(1:NR-1) = initial_condition;
c0(NR) = boundary_condition;
dr = 1/M; %\Delta r
r = 0:dr:1; % dimensionless radius
D = 1.5e-13; %cm^2/s
R = 25e-4; %cm;
alpha0 = D/R^2;

if DCASE == 1
    timevector = [0 0.05 0.1 0.15 0.25 1]./alpha0;
    titlestring = 'Diffusion Case I';
elseif DCASE == 2
    timevector = [0 0.05 0.1175 0.14 0.15 0.2]./alpha0;
    titlestring = 'Diffusion Case II';
elseif DCASE == 3
    timevector = [0 0.05 0.1 0.25 0.5 1]./alpha0;
    titlestring = 'Diffusion Case III';
elseif DCASE == 4
    timevector = [0 0.02 0.05 0.07 0.1 0.28]./alpha0;
    titlestring = 'Diffusion Case IV';
end
    
%% ODE solver (ode45) call
[time,concentration] = ode45(@(t,c)FD_spheres_variable_diffusivity(t,c,DCASE,dr,NR,alpha0),timevector,c0);
c=concentration';
% c = c(r,t) matrix

%% Plotting Concentration Profiles
figure(1) % Plotting c vs. r
hold on
for tt = 1:length(time)
plot(r,c(:,tt),'displayname', ['\alpha_0 t = ' num2str(time(tt)*alpha0)])
end
hold off
xlabel('r')
ylabel('c')
legend('-Dynamiclegend','Location','NorthEast')
title(titlestring)

end
