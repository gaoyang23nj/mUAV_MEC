function Res_Total_Energy = GetTotalEnergy(Num_User, Num_UAV, ck_Rate, Given_TAU_umn,Given_L_un,Given_P_un,Given_Q_mn_x,Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    %% Params
    INIT_PARAMS_K;

    Final_Check = 0;

    %% Main Check Proc
    disp("CheckProc...")
    % % max velocity, 30m/s
    % v_max = 30;
    % % min distance between two UAVs, 50m
    % d_min = 100;

%     %% [0] Get the Target value
%     ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
%     ck_Rate = ck_Rate * Bandwidth;
%     [ck_Target,ck_Delay_Utility,ck_real_Delay_Utility,ck_prop_offload] = GetTargetValue(ck_Rate, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, N, Num_User, Num_UAV)
% 


    %% [1] Energy: total (local computation + communication) energy in User
    disp('[1] Energy: total (local computation + communication) energy in User');
    ck_energy_user_loc_comp = power(Given_L_un,3) * diag(power(Task_Bit_Vec * c_u, 3) * kappa_user /(Delta*Delta));
    ck_energy_user_comm = Given_TAU_umn * (Matrix_Replicate_10_40') .* Given_P_un * Delta;
    ck_energy_user = sum(ck_energy_user_loc_comp + ck_energy_user_comm, 1);
    % Res_Check_energy_user = (ck_energy_user <= E_user_max);
    % [tmp_num_row, tmp_num_col] = size(Res_Check_energy_user);
    % if sum(Res_Check_energy_user(:)) == tmp_num_row * tmp_num_col
    %     fprintf('<%d Pass....%d %d\n',E_user_max, tmp_num_row, tmp_num_col);
    % else
    %     fprintf(' %d', ck_energy_user);
    %     fprintf('\n');
    %     fprintf('<%d No Pass! sum:%d, size:%d %d\n',E_user_max, sum(Res_Check_energy_user(:)), tmp_num_row, tmp_num_col);
    %     Final_Check = Final_Check + 1;
    % end
    

    %% [2] Energy: Propulsion Energy in UAV
    disp('[2] Energy: Propulsion Energy in UAV');
    ck_dist_nplus1_x = [Given_Q_mn_x; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_x];
    ck_dist_nplus1_y = [Given_Q_mn_y; Given_Qinit_mn_y] - [Given_Qinit_mn_y; Given_Q_mn_y];
    ck_velocity_indelta = sqrt(power(ck_dist_nplus1_x, 2) + power(ck_dist_nplus1_y, 2))/Delta;

    ck_prop_energy_term1 = (1 + power(ck_velocity_indelta / prop_param_Utip , 2) * 3) * prop_param_P0;
    ck_prop_energy_term2 = power(ck_velocity_indelta, 3) * prop_param_d0 * prop_param_rho * prop_param_A * prop_param_s * 0.5;
    ck_prop_energy_term3_1 = sqrt(1 + power(ck_velocity_indelta / prop_param_v0, 4)* 0.25);
    ck_prop_energy_term3_2 = power(ck_velocity_indelta / prop_param_v0, 2) * 0.5;
    ck_prop_energy_term3 = sqrt(ck_prop_energy_term3_1 - ck_prop_energy_term3_2) * prop_param_Pi;
    ck_prop_energy = sum(ck_prop_energy_term1 + ck_prop_energy_term2 + ck_prop_energy_term3, 1) * Delta;
    % Res_Check_energy_prop_uav = (ck_prop_energy <= E_uav_prop_max);
    % fprintf(' %d', ck_prop_energy);
    % fprintf('\n');
    % [tmp_num_row, tmp_num_col] = size(Res_Check_energy_prop_uav);
    % if sum(Res_Check_energy_prop_uav(:)) == tmp_num_row * tmp_num_col
    %     fprintf('<%d Pass....%d %d\n', E_uav_prop_max, tmp_num_row, tmp_num_col);
    % else
    %     fprintf('<%d No Pass! sum:%d, size:%d %d\n', E_uav_prop_max, sum(Res_Check_energy_prop_uav(:)), tmp_num_row, tmp_num_col);
    %     Final_Check = Final_Check + 1;
    % end

    %% [3] Energy: Offloading Computation Energy in UAV
    disp('[3] Energy: Offloading Computation Energy in UAV');
    % accurate_Rate
%     ck_uav_x = Given_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
%     ck_uav_y = Given_Q_mn_y * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
%     ck_denomin = power(ck_uav_x, 2) + power(ck_uav_y, 2) + H*H;
%     % Rate hat
%     ck_Rate_hat_lg = ((Given_P_un * rho * Matrix_Replicate_10_40) ./ ck_denomin) * Matrix_Replicate_4_40' + Sigma2;
%     ck_Rate_hat = log2(ck_Rate_hat_lg) .* Bandwidth;
%     % Rate tilde
%     ck_Rate_tilde_lg = ((Given_P_un .* rho * Matrix_Replicate_10_40) ./ ck_denomin) * Matrix_delete_user + Sigma2;
%     ck_Rate_tilde = log2(ck_Rate_tilde_lg) .* Bandwidth;
%     ck_Rate = ck_Rate_hat * Matrix_Replicate_4_40 - ck_Rate_tilde;
    ck_energy_off_comp = power(ck_Rate,3) * kappa_uav * power(c_u, 3) * Delta .* Given_TAU_umn * (Matrix_Replicate_4_40');
    ck_E_uav_comp = sum(ck_energy_off_comp, 1);
    % Res_Check_energy_comp_uav = (ck_E_uav_comp <= E_uav_OE_max);
    % [tmp_num_row, tmp_num_col] = size(Res_Check_energy_comp_uav);
    % if sum(Res_Check_energy_comp_uav(:)) == tmp_num_row * tmp_num_col
    %     fprintf('<%d Pass....%d %d\n',E_uav_OE_max, tmp_num_row, tmp_num_col);
    % else
    %     fprintf(' %d', ck_E_uav_comp);
    %     fprintf('\n');  
    %     fprintf('<%d No Pass! sum:%d, size:%d %d\n',E_uav_OE_max, sum(Res_Check_energy_comp_uav(:)), tmp_num_row, tmp_num_col);
    %     Final_Check = Final_Check + 1;
    % end
    Res_Total_Energy = sum(ck_energy_user) + sum(ck_prop_energy) + sum(ck_E_uav_comp);


end

