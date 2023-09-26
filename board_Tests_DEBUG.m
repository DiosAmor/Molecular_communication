%% Ack
% Special thanks to Birkan
% My codes are based on his work

%% Utilizing method
% Please run the simulator by 'section' in which you want to collect the data

%% Run via Simulator (Poiseuille flow)(for plane source)
clc;
clear;
close all;

% system parameter (micrometer)
% for x, r, phi dimension
% please check the attached image file to check the system model

L=800; R=15; % L: tube length, R: tube radius
 

tx_pts                     = [0 R];     % origin, radius of source plane
tube_pts                   = [0 L R];   % x_min, x_max, radius of tube
receiver_pts               = [L R];   % x, radius of receiving 

            
D_inMicroMeterSqrPerSecond = 100;               % Diffusion coefficient (um^2/s)
average_velocity_inMicroMeterPerSecond = 1000;  % Average velocity (um/s)

alpha = D_inMicroMeterSqrPerSecond*L/(average_velocity_inMicroMeterPerSecond*R^2)

% The number of transmitting molecules =
% delta_sourceplane_radial * delta_sourceplane_angle * ntx_prUnitsource

sim_params.delta_t                      = 5*10^(-4);     %simulation step time (second)
sim_params.delta_sourceplane_radial     = 10;            %source plane radial step 
sim_params.delta_sourceplane_angle      = 100;           %source plane angular step
sim_params.ntx_prUnitsource             = 1;             %transmit molecules per source plane unit
sim_params.tend                         = 2;             %End time of simulation from 0s

replication = 10;

[mol_position_totalstep,mol_arrive_count_avg] = sim_replicator_Poiseuilleflow( ...
   tx_pts,                      ...% Tx points in 3D environment (origin, radius of source plane)
   tube_pts,                    ...% tube size in 3D environment (x_min, x_max, radius of tube)
   receiver_pts,                ...% receiving plane size in 3D environment (x, radius of receiving plane)
   D_inMicroMeterSqrPerSecond,  ...% Diffusion coefficient
   average_velocity_inMicroMeterPerSecond, ...% average velocity of Poiseuille flow
   replication,                 ...% replication number
   sim_params );                   % Simulation parameters

save('data/Poiseuilleflow_V_1000_replication_10_L_800_R_15.mat'); % Change the filename to what you want
fprintf('\n##########################\nRecording Props [done]');

%% Run via Simulator (Uniform flow)(for plane source)
clc;
clear;
close all;

% system parameter (micrometer)
% for x, r, phi dimension
% please check the attached image file to check the system model

L=800; R=15; % L: tube length, R: tube radius

tx_pts                     = [0 R];     % origin, radius of source plane
tube_pts                   = [0 L R];   % x_min, x_max, radius of tube
receiver_pts               = [L L R];   % x_min, x_max, radius of receiving 
            
            
D_inMicroMeterSqrPerSecond = 700;               % Diffusion coefficient (um^2/s)
velocity_inMicroMeterPerSecond = 1000;  % Average velocity (um/s)

alpha = D_inMicroMeterSqrPerSecond*L/(velocity_inMicroMeterPerSecond*R^2)

% The number of transmitting molecules =
% delta_sourceplane_radial * delta_sourceplane_angle * ntx_prUnitsource

sim_params.delta_t                      = 5*10^(-4);     %simulation step time (second)
sim_params.delta_sourceplane_radial     = 10;            %source plane radial step 
sim_params.delta_sourceplane_angle      = 100;           %source plane angular step
sim_params.ntx_prUnitsource             = 1;             %transmit molecules per source plane unit
sim_params.tend                         = 2;             %End time of simulation from 0s

replication = 10;

[mol_position_totalstep,mol_arrive_count_avg] = sim_replicator_Uniformflow( ...
   tx_pts,                      ...% Tx points in 3D environment (origin, radius of source plane)
   tube_pts,                    ...% tube size in 3D environment (x_min, x_max, radius of tube)
   receiver_pts,                ...% receiving plane size in 3D environment (x, radius of receiving plane)
   D_inMicroMeterSqrPerSecond,  ...% Diffusion coefficient
   velocity_inMicroMeterPerSecond, ...% velocity of uniform flow
   replication,                 ...% replication number
   sim_params );                   % Simulation parameters

save('data/Uniformflow_V_1000_replication_10_L_800_R_15.mat'); % Change the filename to what you want
fprintf('\n##########################\nRecording Props [done]');
