function [t,veh_X,veh_V,veh_S] = sim_ring(p,accel_func,veh_n,ring_length,sim_length,dt,print_progress)
% (C) by Benjamin Seibold, Rabie Ramadan, George Gunter

noise = 3e-1; % magnitude of velocity noise in each step
noise0 = 1e-1*0; % magnitude of velocity noise initially
t_noise = 100; % time after which noise is stopped
integrator = 'ballistic';
%integrator = 'Euler';

% Simulation parameters
tf = sim_length; % final time
L = ring_length; % road length

acceleration = @(s,ds,v) accel_func(p,s,ds,v);

% Explanation: Input parameters are:
% s = gap to lead vehicle (x_lead-x_self-length)
% ds = velocity difference from lead vehicle
%      (v_lead-v_self) = -approaching rate [note the difference to IDM]
% v = vehicle velocity
% Output is vehicle acceleration

%% Simulation initialization
%rng(0)
n_steps = ceil(tf/dt); dt = tf/n_steps;
%s_eq = S_eq(v_eq); % corresponding equilibrium gap
%d_eq = s_eq+len; % equilibrium spacing
d_eq = L/veh_n; % road length
s_eq = d_eq;

% Search for velocity eq:

fun = @(v) acceleration(s_eq,0,v);
v_eq = fzero(fun,[0+1e-4,30]); %find equilibrium speed

veh_x = linspace(d_eq,L,veh_n);        % initial vehicle positions
veh_v = veh_x*0+v_eq;                  % initial vehicle velocities
veh_X = veh_x'*ones(1,n_steps+1);      % array storing time-history of positions
veh_V = veh_v'*ones(1,n_steps+1);      % array storing time-history of velocities
veh_S = ones(veh_n,n_steps+1)*s_eq;    % array storing time-history of gaps
veh_v = veh_v+randn(1,veh_n)*noise0;   % add noise to initial velocities

%% Time loop
t = linspace(0,tf,n_steps+1); % time vector
if(print_progress)
    fprintf('_______20%%_______40%%_______60%%_______80%%______100%%\n')
end
for j = 1:n_steps
    if(print_progress)
        if floor(j/n_steps*50)>floor((j-1)/n_steps*50) % progress marker
            fprintf('^')
        end
    end
    % Update step
    veh_s = [veh_x(2:end),veh_x(1)+L]-veh_x;         % vehicle gaps
    veh_ds = veh_v([2:end,1])-veh_v;                     % velocity differences
    dx = veh_v; dv = acceleration(veh_s,veh_ds,veh_v);   % rates of change
    if strcmp(integrator,'ballistic')
        veh_x = veh_x+max(dt*dx+dt^2/2*dv,0);
        veh_v = veh_v+dt*dv;                             % ODE update (ballistic)
    else
        veh_x = veh_x+dt*dx; veh_v = veh_v+dt*dv;        % ODE update (Euler)
    end
    if t(j)<t_noise
        veh_v = veh_v+sqrt(dt)*randn(1,veh_n)*noise;        % add noise to velocity
    end
    veh_v = max(veh_v,0);                                % cap vehicle velocity from below
    
    veh_X(:,j+1) = veh_x';                               % save current state to array
    veh_V(:,j+1) = veh_v';                               % save current state to array
    veh_S(:,j+1) = veh_s';                               % save current state to array
end
fprintf('\n')
