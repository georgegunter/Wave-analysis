%MICRO_MODEL_STABILITY
%   Considers the Aw-Klar-Materne-Rascle model, which is a
%   combination of the follow-the-leader model and the
%   optimal-velocity model:
%      x_j'' = b*(x_{j+1}'-x_j')/(x_{j+1}-x_j)
%            + a*(V(x_{j+1}-x_j)-x_j')
%   The model is initialized with a uniform initial state, plus
%   a small perturbation. Its linearization around this uniform
%   state is found, and the eigenvalues computed. Then, a
%   linear feedback control is applied to the velocity of the
%   first vehicle. The shifted poles of the linearized system
%   are shown, as well as the trajectories of all vehicles for
%   the original system and the feedback-controlled system.
%
%   (C) 2015/11/17 by Benjamin Seibold and Raphael Stern
clear;close all; 
%------------------------------------------------------------------------
%% Model parameters
%------------------------------------------------------------------------
L = 230; % length of circular road
n = 22; % number of vehicles
d0 = 2.23; % reference distance [m]
vm = 9.72; % maximum velocity [m/s]
lv = 4.5; % vehicle length [m]
V = @(d) vm*(tanh(d/d0-2)+tanh(2))/(1+tanh(2)); % optimal velocity function
b = 20; % follow-the-leader strength [m^2/s]
a = 0.5; % optimal velocity strength [1/s]
%b = 0.5; % follow-the-leader strength [m^2/s]
%a = 2.2; % optimal velocity strength [1/s] 
tf = 300; % final time of trajectory computation
sigma = .25; % noise level [m/s^2]
dt_noise = 1; % time between noise events [s]
dt = .1; % time step for output [s]
B = (1:2*n)'==2; % control matrix (can affect only velocity of vehicle 1)
t_av_on = 30; % time at which the AV turns on the feedback control

% NOTE: If you want the vehicle to have control increase the values of
% c_control, which is currently set to 0.
%------------------------------------------------------------------------
%% Derived parameters and functions
%------------------------------------------------------------------------
d = L/n-lv; % initial distance of adjacent vehicles
dV = @(d) (V(d+1e-7)-V(d-1e-7))/2e-7;
fprintf('Microscopic model (AKMR) on circular road of length %g.\n',L)
fprintf('%d vehicles; vehicle spacing: %0.3g\n',n,d)
fprintf('opt. vel.: V(d) = %0.3g; V''(d) = %0.3g\n',V(d),dV(d))
dist = @(x) [x(3:2:end);x(1)+L]-x(1:2:end)-lv; % distance between vehicles
% ODE system right hand side function
R = @(~,x) reshape([x(2:2:end),b*(x([4:2:end,2])-x(2:2:end))./dist(x).^2 + ...
    a*(V(dist(x))-x(2:2:end))]',[],1);
R_c = @(~,x) R(0,x)+B*(F*(x-V(d))); % feedback-controlled right hand side

%------------------------------------------------------------------------
%% Set up linearized system matrices and calculate eigenmodes
%------------------------------------------------------------------------
% Solution x_j(t) = dj+ut + y_j(t)
% x_j'' = b*(x_{j+1}'-x_j')/(x_{j+1}-x_j-lv)^2 + a*(V(x_{j+1}-x_j-lv)-x_j')
% x' = u+y'   and   V(d-lv) = u
% y_j'' = b*(y_{j+1}'-y_j')/(d+y_{j+1}-y_j-lv)^2 + a*(V(d+y_{j+1}-y_j-lv)-u-y_j')
% y_j'' = b/(d-lv)^2*(y_{j+1}'-y_j') + a*(V'(d-lv)*(y_{j+1}-y_j)-y_j')
% y_j'' = b/(d-lv)^2*y_{j+1}' - (a+b/(d-lv)^2)*y_j' + a*V'(d-lv)*(y_{j+1}-y_j)

cb = b/d^2; cc = a*dV(d); ca = a+cb; % model parameters
av = (1+(-1).^(1:2*n))/2; % alternating vector
A = diag(-ca*av)+diag(-cc*av(2:end),-1)+...       % linearized
    diag(cc+(1-cc)*av(2:end),1)+diag(cb*av(3:end),2); % system
A(end,1:2) = [cc cb];                      % matrix (circular)

%full(A)

%[EV,D] = eig(A); lambda = diag(D);
lambda = eig(A); % eigenvalues of linearized system
ind_unstable = find(real(lambda)>1e-14); % indices of unstable modes
fprintf('number of unstable modes: %d\n',numel(ind_unstable))


c_control = 10.0; % strength of feedback control, gain
F = B'*0; F(2) = F(2)-c_control; % feedback matrix
A_c = A+B*F; % feedback-controlled matrix of linearized system
lambda_c = eig(A_c); % eigenvalues of feedback-controlled matrix
fprintf('feedback-controlled system: unstable modes: %d\n',...
    sum(real(lambda_c)>1e-14))

%------------------------------------------------------------------------
%% Run computation of microscopic model (original; not linearized)
%------------------------------------------------------------------------
p = 0:(d+lv):(d+lv)*(n-1); % initial positions of vehicles
v = p*0+V(d); % initial velocities of vehicles
%v = v+(rand(size(v))*2-1)*d*.1; % add noise to initial velocities
x0 = reshape([p;v],[],1); % initial state vector
% Time stepping
dt = dt_noise/ceil(dt_noise/dt); % integer number of steps per subinterval
t_step = 0:dt:dt_noise; % time steps in each subinterval
n_sub = ceil(tf/dt_noise); % number of sub-intervals
sigma_dt = sqrt(dt_noise)*sigma; % noise for sub-interval

rng(1) % seed random number generator
t = 0; Y = x0'; Y_c = Y; % initial time vector and state matrix
for j = 1:n_sub % loop over sub-intervals
    % add noise
    %%%%%%%%%%%%%%%%%%%%%%%
    %Y(end,2:2:end) = Y(end,2:2:end)+randntrunc(1,n,3)*sigma_dt;
    %Y_c(end,2:2:end) = Y_c(end,2:2:end)+[0,randntrunc(1,n-1,3)]*sigma_dt;
    %%%%%%%%%%%%%%%%%%%%%%%

    noise = randntrunc(1,n,3)*sigma_dt; % compute the noise to be added in each sub interval
    noise_c = noise; noise_c(1)=0; % re-use this noise but do not apply noise to the first vehicle (AV)
    Y(end,2:2:end) = Y(end,2:2:end)+noise; % add noise to positions in uncontrolled case
    Y_c(end,2:2:end) = Y_c(end,2:2:end)+noise_c; % add noise to positions in controlled case
    
    % advance ODEs
    t_sub = (j-1)*dt_noise+t_step; % sub-interval time steps
    [~,Y_sub] = ode45(R,t_sub,Y(end,:)); % compute vehicle trajectories
    options = odeset('Events', @endEvent); % set the options with the myEvent function
    if t_sub<t_av_on % before the AV is switched on, use the un-controlled evolution
        [~,Y_c_sub] = ode45(R,t_sub,Y_c(end,:), options);
    else
        [~,Y_c_sub] = ode45(R_c,t_sub,Y_c(end,:), options); % compute veh. trajectories
    end
    
    t = [t,t_sub(2:end)]; % append states
    Y = [Y;Y_sub(2:end,:)]; Y_c = [Y_c;Y_c_sub(2:end,:)];
end
P = mod([Y(:,1:2:end),Y(:,1:2:end)+L,Y(:,1:2:end)+2*L],3*L)-L;
P(P<-.5*L|P>1.5*L) = nan;
P_c = mod([Y_c(:,1:2:end),Y_c(:,1:2:end)+L,Y_c(:,1:2:end)+2*L],3*L)-L;
P_c(P_c<-.5*L|P_c>1.5*L) = nan;

loc = Y(:,1:2:end); % location of each car at every point in time
loc_c = Y_c(:,1:2:end); % location of each car when control is on
headway_av = mod(loc(:,2)-loc(:,1)+L, L); % headway of AV
headway_av_c = mod(loc_c(:,2)-loc_c(:,1)+L, L); % headway of AV when control is on


% compute headways of each vehicle (except for the AV)
h = []; % the headway matrix
h_c = []; % headway of controlled system
for i=1:n % loop over each vehicle to compute the headway at each instance in time
h(:,i) = mod(loc(:,max(1, mod(i+1, n+1)))-loc(:,i)+L, L);
h_c(:,i) = mod(loc_c(:,max(1, mod(i+1, n+1)))-loc_c(:,i)+L, L);
end
h(:,1) = []; % delete the AV
h_c(:,1) = [];
min_h = min(h, [], 2);
min_h_c = min(h_c, [], 2);
max_h = max(h, [], 2);
max_h_c = max(h, [], 2);


%------------------------------------------------------------------------
%% Plot results
%------------------------------------------------------------------------
% clf
% subplot(3,2,1)
% plot(t,P,'b-',t,P(:,1:n:end),'r-')
% axis([0 tf 0 L])
% xlabel('time t'), ylabel('position on road x')
% title('vehicle trajectories of original system')
% 
% subplot(3,2,2)
% plot(t,P_c,'b-',t,P_c(:,1:n:end),'r-')
% axis([0 tf 0 L])
% xlabel('time t'), ylabel('position on road x')
% title('vehicle trajectories of feedback-controlled system')
% 
% subplot(3,2,3)
% plot(t, headway_av, 'r-', [0, tf], [L/n, L/n], 'k--', t, min_h, 'b-', t, max_h, 'b-')
% axis([0, tf, 0, max(max(headway_av), max(headway_av_c))*1.2])
% xlabel('time t'), ylabel('headway (m)')
% title('headway of vehicles')
% 
% subplot(3,2,4)
% plot(t, headway_av_c, 'r-', [0, tf], [L/n, L/n], 'k--', t, min_h_c, 'b-', t, max_h, 'b-')
% axis([0, tf, 0, max(max(headway_av), max(headway_av_c))*1.2])
% xlabel('time t'), ylabel('headway (m)')
% title('headway of vehicles')
% 
% subplot(3,2,5)
% p = sort(P(end,:)); p = p(~isnan(p)); % positions at final time
% p_mid = (p(1:end-1)+p(2:end))/2; % mid points between vehicles
% rho = 1./diff(p); % vehicle density
% p_c = sort(P_c(end,:)); % controlled system
% veh1 = p_c(1:n:end); veh1 = veh1(0<=veh1&veh1<L); % pos. of vehicle 1
% p_c = p_c(~isnan(p_c)); % controlled system: positions at final time
% p_mid_c = (p_c(1:end-1)+p_c(2:end))/2; % mid points between vehicles
% rho_c = 1./diff(p_c); % vehicle density
% ax = [0 L 0 max([rho rho_c])*1.05];
% plot(p_mid,rho,'k-',p_mid_c,rho_c,'r-',veh1*[1 1],ax(3:4),'r:')
% axis(ax)
% legend('original system','controlled system')
% xlabel('position on road x'), ylabel('vehicle density \rho')
% title(sprintf('density at final time (t=%g)',tf))
% 
% subplot(3,2,6)
% lar = real(lambda); lai = imag(lambda);
% ax = [min(lar),max(lar),min(lai),max(lai)];
% ax = ax+[[-1 1]*(ax(2)-ax(1)),[-1 1]*(ax(4)-ax(3))]*.05;
% plot([0 0],ax(3:4),'k-',...
%     lar,lai,'b.',lar(ind_unstable),lai(ind_unstable),'bo',...
%     real(lambda_c),imag(lambda_c),'r.')
% axis(ax)
% xlabel('real part'), ylabel('imaginary part')
% title('eigenvalues of linearized system in complex plane')




clf
subplot(2,1,1)
plot(t,P,'b-',t,P(:,1:n:end),'r-')
axis([0 tf 0 L])
xlabel('Time t (s)'), ylabel('Position on road x (m)')
title('Vehicle trajectories of original system')

subplot(2,1,2)
plot(t,P_c,'b-',t,P_c(:,1:n:end),'r-')
axis([0 tf 0 L])
xlabel('Time t (s)'), ylabel('Position on road x (m)')
title('Vehicle trajectories of feedback-controlled system')

%subplot(2,2,3)
%p = sort(P(end,:)); p = p(~isnan(p)); % positions at final time
%p_mid = (p(1:end-1)+p(2:end))/2; % mid points between vehicles
%rho = 1./diff(p); % vehicle density
%p_c = sort(P_c(end,:)); % controlled system
%veh1 = p_c(1:n:end); veh1 = veh1(0<=veh1&veh1<L); % pos. of vehicle 1
%p_c = p_c(~isnan(p_c)); % controlled system: positions at final time
%p_mid_c = (p_c(1:end-1)+p_c(2:end))/2; % mid points between vehicles
%rho_c = 1./diff(p_c); % vehicle density
%ax = [0 L 0 max([rho rho_c])*1.05];
%plot(p_mid,rho,'k-',p_mid_c,rho_c,'r-',veh1*[1 1],ax(3:4),'r:')
%axis(ax)
%legend('original system','controlled system')
%xlabel('position on road x'), ylabel('vehicle density \rho')
%title(sprintf('density at final time (t=%g)',tf))

%subplot(2,2,4)
%lar = real(lambda); lai = imag(lambda);
%ax = [min(lar),max(lar),min(lai),max(lai)];
%ax = ax+[[-1 1]*(ax(2)-ax(1)),[-1 1]*(ax(4)-ax(3))]*.05;
%plot([0 0],ax(3:4),'k-',...
%    lar,lai,'b.',lar(ind_unstable),lai(ind_unstable),'bo',...
%    real(lambda_c),imag(lambda_c),'r.')
%axis(ax)
%xlabel('real part'), ylabel('imaginary part')
%title('eigenvalues of linearized system in complex plane')

%------------------------------------------------------------------------
%% Save results
%------------------------------------------------------------------------
%save('trajectories','t','Y','Y_c')





