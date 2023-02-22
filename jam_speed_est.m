function [jam_speed_est] = jam_speed_est(p,accel_func,veh_n,ring_length,sim_length,dt,print_progress)
%% Run a sim:
[t,veh_X,veh_V,~] = ...
    sim_ring(p,accel_func,veh_n,ring_length,sim_length,dt,print_progress);
%% Get agg:
wave_start_point = floor(sim_length*dt)/2;%Look at last half
t = t(wave_start_point:end);
veh_X = veh_X(:,wave_start_point:end);
veh_V = veh_V(:,wave_start_point:end);


[~,~,instances] = ...
    extract_agg_data(t,veh_X,veh_V,ring_length/2,30,ring_length);

instances = sortrows(instances,1);
%% estimate jam speed:
counts = (diff(instances(:,1))).^-1; %sec/veh
flows = counts*3600; %veh/hr
speeds_agg = instances(1:end-1,2); %m/s
densities_agg = flows./(speeds_agg*3600); %veh/m

jam_speed_est = polyfit(densities_agg*1000,flows,1);
jam_speed_est = jam_speed_est(1);

end

