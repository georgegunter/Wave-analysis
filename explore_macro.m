clear;
clc;
%% Instantiate sim params: IDM
accel_func = @(p,s,ds,v) p(1)*( 1 - (v/p(3)).^p(4) - ...
    ((p(5) + max(p(6)*v + v.*(-ds)/(2*sqrt(p(1)*p(2))),0))./s).^2 );
p = [1.3,2.0,30,4,2.0,1.0];

%% OVM-FTL:
accel_func = @(p,s,ds,v)  ...
    p(1)*(p(3)/2*(1 - cos(pi*(s-p(4))/(p(5) - p(4))))-v) + p(2)*(ds);

p = [.1,5,9.8,2,15];

%% Bando OVM-FTL:
accel_func = @(p,s,ds,v)  ...
    p(1)*(p(3)*(tanh(s./p(4)-p(5))+tanh(p(5)))/(1+tanh(p(5)))-v) + p(2)*((ds)./(s.^2));

p = [.5,20,30.0,12,2];

%% Sim params:
veh_n = 21;
ring_length = 154.14; %Ared exp G (without vehicle lengths)
sim_length = 900;
dt = 1/30;
print_progress = true;
%% Run sim:

[times,veh_X,veh_V,veh_S] = sim_ring(p,accel_func,veh_n,ring_length,sim_length,dt,print_progress);

%% Select relevant range:
n=12000;
times = times(n:end);
veh_X = veh_X(:,n:end);
veh_V = veh_V(:,n:end);
veh_S = veh_S(:,n:end);

%% Extract Macro estimates:

[MRho,MQ,jam_speed,mt,mx] = ...
    get_macro_filter_data_ring(times,ring_length,veh_X,veh_V);
[MT,MX] = meshgrid(mt,mx);
MU = MQ./MRho;

disp('Finished extracting macro quantities.')
%% Look at numerical gradients:
dt = mt(2) - mt(1);
dx = mx(2) - mx(1);

[dU_dt,dU_dx] = gradient(MU,dt,dx);
[dQ_dt,dQ_dx] = gradient(MQ,dt,dx);
[dRho_dt,dRho_dx] = gradient(MRho,dt,dx);

% [dQ_dx,dQ_dt] = gradient(MQ,dx,dt);
% [dRho_dx,dRho_dt] = gradient(MRho,dx,dt);

%% plot continuity:

cont = dRho_dt + dQ_dx; %Should be zero?
figure()
contourf(cont)
colorbar()


%% Pre-processing:
num_vehicles = length(veh_X(:,1));
num_samples = length(veh_X(1,:));

veh_X_ring = mod(veh_X,ring_length);

veh_T = ones(num_vehicles,num_samples);
for i = 1:num_vehicles
    veh_T(i,:) = times;
end


veh_T_reshaped = reshape(veh_T,num_vehicles*num_samples,1);
veh_X_reshaped = reshape(veh_X_ring,num_vehicles*num_samples,1);
veh_V_reshaped = reshape(veh_V,num_vehicles*num_samples,1);

%% Plot the time-space diagram vs. macro speed estimate:
figure()
hold on
scatter(veh_T_reshaped,mod(veh_X_reshaped,ring_length),1.0,veh_V_reshaped)
xlim([times(1),times(end)]);
ylim([0,ring_length]);
box on
ylabel('Position [m]')
xlabel('Time [s]')
set(gca,'FontSize',24)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',5)
title('Space-Time')
colorbar()
%%
figure()
contourf(MT,MX,MU)
ylabel('Position [m]')
xlabel('Time [s]')
set(gca,'FontSize',24)
set(gca,'FontWeight','bold')
set(gca,'LineWidth',5)
title('Macroscopic Speed Estimate')
colorbar()

