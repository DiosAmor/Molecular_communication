clc;
clear;
close all;

load('data/Poiseuilleflow_V_1000_replication_10_L_800_R_15.mat');
micro_plot = figure();
set(micro_plot, 'Units', 'centimeters')
set(micro_plot, 'Position', [0 0 25 25]);

index=1;
for t=1:200:2000
    subplot(5,2,index);
    x = mol_position_totalstep(:,1,t);
    y = mol_position_totalstep(:,2,t);
    scatter(x,y);
    xlim([0 L]);
    ylim([-R R]);
    xlabel('x axis (um)');ylabel('y axis (um)');
    a = 'Poiseuille flow: ';
    b =  ['molecules postions at time ' num2str((t-1)*sim_params.delta_t) 's'];
    title({a,b});
    index=index+1;
end


load('data/Uniformflow_V_1000_replication_10_L_800_R_15.mat');
mol_arrive_count_avg_micro = mol_arrive_count_avg;
macro_plot = figure();
set(macro_plot, 'Units', 'centimeters')
set(macro_plot, 'Position', [0 0 25 25]);

index=1;
for t=1:200:2000
    subplot(5,2,index);
    x = mol_position_totalstep(:,1,t);
    y = mol_position_totalstep(:,2,t);
    scatter(x,y);
    xlim([0 L]);ylim([-R R]);
    xlabel('x axis (um)');ylabel('y axis (um)');
    a = 'Uniform flow: ';
    b =  ['molecules postions at time ' num2str((t-1)*sim_params.delta_t) 's'];
    title({a,b});
    index=index+1;
end