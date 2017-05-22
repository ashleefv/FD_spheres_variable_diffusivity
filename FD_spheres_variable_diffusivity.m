function dcdt = FD_spheres_variable_diffusivity(t,c,DCASE,dr,NR,alpha0)
%   Finite difference discritization scheme for spheres
%   with variable diffusivity as described in the reference (scheme 1).
%
%% --- Reference --- 
%  A. N. Ford Versypt & R. D. Braatz, Analysis of finite difference 
%  discretization schemes for diffusion in spheres with variable 
%  diffusivity, Computers & Chemical Engineering, 71 (2014) 241-252, 
%  https://doi.org/10.1016/j.compchemeng.2014.05.022
% 
%% --- Input ---
% 	alpha(NR)   D(r,t)/R^2 for the species
% 	dr          Spatial discretization size
% 	c(NR)       ODE solution vector
%	NR          Number of spatial discretizations
%% --- Output ---
%   dcdt(NR)      ODE derivative vector updated with diffusion contribution

%% Nested function defining the diffusivity for each specified DCASE
    function alpha = diffusivity(DCASE,t,NR,dr,c,alpha0)
        % Generates alpha(r) at each time t
        alpha = zeros(1,NR); % initialize an empty vector for alpha
        
        if DCASE == 1
            % Case I: constant diffusivity with alpha = alpha0
            alpha = alpha0*ones(1,NR);
        elseif DCASE == 2
            % Case II: time-dependent diffusivity
            tau = 0.1175/alpha0;
            k = 10/tau; 
            if t < tau
                alpha = alpha0*ones(1,NR); 
            else
                alpha = alpha0*exp(k*(t - tau))*ones(1,NR); 
            end
        elseif DCASE == 3
            % Case III: spatially dependent diffusivity
            r = 0:dr:1;
            alpha = alpha0*r.^2;
        elseif DCASE == 4
            % CASE IV: concentration-dependent diffusivity
            k = 1;
            alpha = alpha0*exp(k-k*c);
        end
    end

%% Compute alpha for 1:NR for the current time and concentration
alpha = diffusivity(DCASE,t,NR,dr,c,alpha0);

%% Compute the diffusion term at r=0
    for rr = 1
     	dcdt(rr) = 6*alpha(rr)*(c(rr+1)-c(rr))/dr^2;
    end
%% Compute the diffusion terms at 0<r<1
% --- Loop over all interior spatial discretizations 0<r<1 ---
	for rr = 2:NR-1
        i = rr-1; % i indexing starts at 0 and rr indexing in MATLAB starts at 1
%       scheme 1 from the reference
        dcdt(rr) = alpha(rr)/i/2/dr^2*...
        ((i+2)*c(rr+1)-2*i*c(rr)+(i-2)*c(rr-1))...
        +1/2/dr^2*alpha(rr+1)*(c(rr+1)-c(rr))...
        +1/2/dr^2*alpha(rr-1)*(c(rr-1)-c(rr));
    end    
%% Constant boundary condition at r = 1    
    for rr = NR
        dcdt(rr) = 0;
    end
    
    dcdt = dcdt';
end
