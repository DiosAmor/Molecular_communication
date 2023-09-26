function [mol_position_totalstep, mol_arrive_count_avg ] = sim_replicator_Poiseuilleflow( ...
   tx_pts,                      ...% Tx points in 3D environment (origin, radius of source plane)
   tube_pts,                    ...% tube size in 3D environment (x_min, x_max, radius of tube)
   receiver_pts,                ...% receiving plane size in 3D environment (x, radius of receiving plane)
   D_inMicroMeterSqrPerSecond,  ...% Diffusion coefficient
   average_velocity_inMicroMeterPerSecond, ...% average velocity of Poiseuille flow
   replication, ...% replication number
   sim_params )   % Simulation parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates convection & diffusion channel in 3D env
% with a planar source(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% sim_params.delta_t                      = simulation step time;
% sim_params.delta_sourceplane_radial     = source plane radial step ;
% sim_params.delta_sourceplane_angle      = source plane angular step;
% sim_params.ntx_prUnitsource             = transmit molecules per source plane unit;
% sim_params.tend                         = End time of simulation from 0s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

tStart = tic;
[mol_position_totalstep, mol_arrive_count] = sim_3d_PoS_Poiseuilleflow( ...
   tx_pts,                      ...% Tx points in 3D environment (origin, radius of source plane)
   tube_pts,                    ...% tube size in 3D environment (x_min, x_max, radius of tube)
   receiver_pts,                ...% receiving plane size in 3D environment (x, radius of receiving plane)
   D_inMicroMeterSqrPerSecond,  ...% Diffusion coefficient
   average_velocity_inMicroMeterPerSecond, ...% average velocity of Poiseuille flow
   sim_params ); % Simulation parameters

mol_arrive_count_avg = mol_arrive_count/replication;
tElapsed = toc(tStart);
fprintf('\n ###########  replication=%d / %d [DONE] in (%f sec)',1,replication,tElapsed);

for replication_idx=2:replication
    tStart = tic;
    [mol_position_totalstep, mol_arrive_count] = sim_3d_PoS_Poiseuilleflow( ...
       tx_pts,                      ...% Tx points in 3D environment (origin, radius of source plane)
       tube_pts,                    ...% tube size in 3D environment (x_min, x_max, radius of tube)
       receiver_pts,                ...% receiving plane size in 3D environment (x, radius of receiving plane)
       D_inMicroMeterSqrPerSecond,  ...% Diffusion coefficient
       average_velocity_inMicroMeterPerSecond, ...% average velocity of Poiseuille flow
       sim_params ); % Simulation parameters
    mol_arrive_count_avg = mol_arrive_count_avg + mol_arrive_count/replication;
    tElapsed = toc(tStart);
    fprintf('\n ###########  replication=%d / %d [DONE] in (%f sec)',replication_idx,replication,tElapsed);
end

end
