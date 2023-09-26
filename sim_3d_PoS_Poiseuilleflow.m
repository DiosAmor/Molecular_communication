function [mol_position_totalstep, mol_arrive_count] = sim_3d_PoS_Poiseuilleflow( ...
   tx_pts,                      ...% Tx points in 3D environment (origin, radius of source plane)
   tube_pts,                    ...% tube size in 3D environment (x_min, x_max, radius of tube)
   receiver_pts,                ...% receiving plane size in 3D environment (x, radius of receiving plane)
   D_inMicroMeterSqrPerSecond,  ...% Diffusion coefficient
   average_velocity_inMicroMeterPerSecond, ...% average velocity of Poiseuille flow
   sim_params ) % Simulation parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates diffusion channel in 3D env
% with a planar source(s) and an absorbing receiving plane
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


D                       = D_inMicroMeterSqrPerSecond;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = sim_params.tend;
ntx_radial              = sim_params.delta_sourceplane_radial;
ntx_angle               = sim_params.delta_sourceplane_angle;
ntx_unit                = sim_params.ntx_prUnitsource;
V                       = average_velocity_inMicroMeterPerSecond;

%
tube_x_low  = tube_pts(1);
tube_x_high = tube_pts(2);
tube_radius = tube_pts(3);

receiver_x  = receiver_pts(1);
receiver_radius = receiver_pts(2);


% First find the number of simulation steps (time)
sim_step_cnt = round(sim_time_inSeconds / delta_t);

% Standard deviation of step size of movement N(0,sigma)
sigma = (2*D*delta_t)^0.5;

% Rx timeline Records the number of molecules at Receiving plane at each time step 
% First, Distribute molecules uniformly at the source plane
ntx = ntx_radial*ntx_angle*ntx_unit;
mol_position1 = repmat([0 0 0], ntx, 1);
mol_index=1;
interval_radius = tx_pts(2)/ntx_radial;
interval_angle = 2*pi/ntx_angle;

for ir=0:interval_radius:tx_pts(2)-tx_pts(2)/ntx_radial
    for ia=0:interval_angle:2*pi-2*pi/ntx_angle
        for iu=1:ntx_unit
            mol_position1(mol_index,:) = [tx_pts(1) ir*cos(ia) ir*sin(ia)];
            mol_index = mol_index+1;
        end
    end
end
mol_position_totalstep = zeros(ntx,3,sim_step_cnt+1);   % Record traces of molecules by all time
mol_arrive_count = zeros(sim_step_cnt,1);               % Record the number of arrived molecules
mol_position_totalstep(:,:,1) = mol_position1;

% Seocond, Move molecules
for t=2:sim_step_cnt
    % Propagate the molecules via diffusion & drift
    % Poiseuilleflow shows the velocity profile which is fastest at the middle and 0 at the surface.
    mol_displace = normrnd(0, sigma, ntx , 3); % diffusion
    mol_displace_v = [2*V*(1-(mol_position1(:,2).^2+mol_position1(:,3).^2)/tube_radius^2)*delta_t zeros(ntx,2)]; % drift
    
    % check already arrive (Already arrived molecules doesn't move because it is absorbed)
    arrive_mask_already = (mol_position1(:,1) > receiver_x) | mol_position1(:,1) == receiver_x;
    mol_displace(arrive_mask_already,:) = zeros(sum(arrive_mask_already),3);
    mol_displace_v(arrive_mask_already,:) = zeros(sum(arrive_mask_already),3);
    mol_position2 =  mol_position1 + mol_displace + mol_displace_v;

    % Outside the tube 
    outside_tube_x_mask = mol_position2(:,1) < tube_x_low;
    outside_tube_radius_mask = (mol_position2(:,2).^2 + mol_position2(:,3).^2) > tube_radius^2;
    % Roll back the molecules inside the tube
    mol_position2(outside_tube_x_mask | outside_tube_radius_mask, :) = mol_position1(outside_tube_x_mask | outside_tube_radius_mask, :);
    
    % arrive molecules inside receiver
    arrive_mask = (mol_position2(:,1) > receiver_x) & (mol_position2(:,2).^2+mol_position2(:,3).^2 < receiver_radius^2) & (mol_position1(:,1) < receiver_x);
    
    % count arrive molecules
    mol_arrive_count(t) = sum(arrive_mask~=0);
    mol_position1 = mol_position2;
    mol_position_totalstep(:,:,t) = mol_position1;
end



