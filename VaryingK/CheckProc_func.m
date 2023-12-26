function [Final_Check, ck_prop_energy] = CheckProc_func(Num_User, Num_UAV, ck_Rate, Given_TAU_umn,Given_L_un,Given_P_un,Given_Q_mn_x,Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec)
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
    ck_energy_user_loc_comp = power(Given_L_un,3) * diag(power(Task_Bit_Vec * c_u, 3) * kappa_user /(delta*delta));
    ck_energy_user_comm = Given_TAU_umn * (Matrix_Replicate_10_40') .* Given_P_un * delta;
    ck_energy_user = sum(ck_energy_user_loc_comp + ck_energy_user_comm, 1);
    Res_Check_energy_user = (ck_energy_user <= E_user_max);
    [tmp_num_row, tmp_num_col] = size(Res_Check_energy_user);
    if sum(Res_Check_energy_user(:)) == tmp_num_row * tmp_num_col
        fprintf('<%d Pass....%d %d\n',E_user_max, tmp_num_row, tmp_num_col);
    else
        fprintf(' %d', ck_energy_user);
        fprintf('\n');
        fprintf('<%d No Pass! sum:%d, size:%d %d\n',E_user_max, sum(Res_Check_energy_user(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

    %% [2] Energy: Propulsion Energy in UAV
    disp('[2] Energy: Propulsion Energy in UAV');
    ck_dist_nplus1_x = [Given_Q_mn_x; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_x];
    ck_dist_nplus1_y = [Given_Q_mn_y; Given_Qinit_mn_y] - [Given_Qinit_mn_y; Given_Q_mn_y];
    ck_velocity_indelta = sqrt(power(ck_dist_nplus1_x, 2) + power(ck_dist_nplus1_y, 2))/delta;

    ck_prop_energy_term1 = (1 + power(ck_velocity_indelta / prop_param_Utip , 2) * 3) * prop_param_P0;
    ck_prop_energy_term2 = power(ck_velocity_indelta, 3) * prop_param_d0 * prop_param_rho * prop_param_A * prop_param_s * 0.5;
    ck_prop_energy_term3_1 = sqrt(1 + power(ck_velocity_indelta / prop_param_v0, 4)* 0.25);
    ck_prop_energy_term3_2 = power(ck_velocity_indelta / prop_param_v0, 2) * 0.5;
    ck_prop_energy_term3 = sqrt(ck_prop_energy_term3_1 - ck_prop_energy_term3_2) * prop_param_Pi;
    ck_prop_energy = sum(ck_prop_energy_term1 + ck_prop_energy_term2 + ck_prop_energy_term3, 1) * delta;
    Res_Check_energy_prop_uav = (ck_prop_energy <= E_uav_prop_max);
    fprintf(' %d', ck_prop_energy);
    fprintf('\n');
    [tmp_num_row, tmp_num_col] = size(Res_Check_energy_prop_uav);
    if sum(Res_Check_energy_prop_uav(:)) == tmp_num_row * tmp_num_col
        fprintf('<%d Pass....%d %d\n', E_uav_prop_max, tmp_num_row, tmp_num_col);
    else
        fprintf('<%d No Pass! sum:%d, size:%d %d\n', E_uav_prop_max, sum(Res_Check_energy_prop_uav(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

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
    ck_energy_off_comp = power(ck_Rate,3) * kappa_uav * power(c_u, 3) * delta .* Given_TAU_umn * (Matrix_Replicate_4_40');
    ck_E_uav_comp = sum(ck_energy_off_comp, 1);
    Res_Check_energy_comp_uav = (ck_E_uav_comp <= E_uav_OE_max);
    [tmp_num_row, tmp_num_col] = size(Res_Check_energy_comp_uav);
    if sum(Res_Check_energy_comp_uav(:)) == tmp_num_row * tmp_num_col
        fprintf('<%d Pass....%d %d\n',E_uav_OE_max, tmp_num_row, tmp_num_col);
    else
        fprintf(' %d', ck_E_uav_comp);
        fprintf('\n');  
        fprintf('<%d No Pass! sum:%d, size:%d %d\n',E_uav_OE_max, sum(Res_Check_energy_comp_uav(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

    %% [4] Rate: Frequency Needed in UAV
    disp('[4] Rate: Maximum Frequency Needed in UAV');
    freq_needed = ck_Rate * c_u;
    Res_Check_freq = (freq_needed <= CPUFreq_UAV);
    [tmp_num_row, tmp_num_col] = size(Res_Check_freq);
    if sum(Res_Check_freq(:)) == tmp_num_row * tmp_num_col
        fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
    else
        fprintf('sum:%d, size:%d %d\n',sum(Res_Check_freq(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

    %% [5] Mobility: Maximum Velocity of UAV
    disp('[5] Mobility: Maximum Velocity of UAV');
    ck_dist_nplus1_x = [Given_Q_mn_x; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_x];
    ck_dist_nplus1_y = [Given_Q_mn_y; Given_Qinit_mn_y] - [Given_Qinit_mn_y; Given_Q_mn_y];
    % dist_nplus1_x = [Given_Q_mn_x_r; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_x_r];
    % dist_nplus1_y = [Given_Q_mn_y_r; Given_Qinit_mn_y] - [Given_Qinit_mn_y; Given_Q_mn_y_r];
    ck_rel_velocity_indelta = (power(ck_dist_nplus1_x, 2) + power(ck_dist_nplus1_y, 2))/(power(delta*v_max, 2));
    Res_Check_velocity = (ck_rel_velocity_indelta <= 1);
    [tmp_num_row, tmp_num_col] = size(Res_Check_velocity);
    if sum(Res_Check_velocity(:)) == tmp_num_row * tmp_num_col
        fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
    else
        fprintf('sum:%d, size:%d %d\n',sum(Res_Check_velocity(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

    %% [6] Mobility: Minimum Distance between 2 UAVs (safety)
    disp('[6] Mobility: Minimum Distance between 2 UAVs (safety)');
    Matrix_2UAV = zeros(Num_UAV, Num_UAV * 2);
    for i=1:Num_UAV
        for j=1:Num_UAV-1
            tmp = zeros(Num_UAV, 1);
            tmp(i) = 1;
            if j<i
                tmp(j,1) = -1;
            else
                tmp(j+1,1) = -1;
            end
            Matrix_2UAV(:, (i-1)*(Num_UAV-1)+j) = tmp;
        end
    end
    % St_2UAV_term1 = -(power(Given_Q_mn_x_r * Matrix_2UAV, 2) + power(Given_Q_mn_y_r * Matrix_2UAV, 2));
    % % St_2UAV_term2 = ((Given_Q_mn_x_r * Matrix_2UAV) .* (Var_Q_mn_x * Matrix_2UAV)) + ((Given_Q_mn_y_r * Matrix_2UAV) .* (Var_Q_mn_y * Matrix_2UAV));
    % St_2UAV_term2 = ((Given_Q_mn_x_r * Matrix_2UAV) .* (Given_Q_mn_x_r * Matrix_2UAV)) + ((Given_Q_mn_y_r * Matrix_2UAV) .* (Given_Q_mn_y_r * Matrix_2UAV));
    % St_2UAV = (St_2UAV_term1 + St_2UAV_term2 * 2)/1e5;
    ck_distance_2UAVs = power(Given_Q_mn_x * Matrix_2UAV, 2) + power(Given_Q_mn_y * Matrix_2UAV, 2);
    Res_Check_Dis2UAVs = (ck_distance_2UAVs >= (d_min*d_min));
    [tmp_num_row, tmp_num_col] = size(Res_Check_Dis2UAVs);
    if sum(Res_Check_Dis2UAVs(:)) == tmp_num_row * tmp_num_col
        fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
    else
        fprintf('sum:%d, size:%d %d\n',sum(Res_Check_Dis2UAVs(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

    %% [7] Bits: Task should be completed
    disp('[7] Bits: Task should be completed');
    ck_offload_bits = ck_Rate .* Given_TAU_umn * delta * (Matrix_Replicate_10_40');
    ck_local_bits = Given_L_un * diag(Task_Bit_Vec);
    ck_Computed_bits = sum(ck_offload_bits + ck_local_bits, 1);
    Res_Check_TaskCompleted = (ck_Computed_bits >= Task_Bit_Vec);
    % [tmp_num_row, tmp_num_col] = size(Res_Check_TaskCompleted);
    if sum(Res_Check_TaskCompleted(:)) == length(Res_Check_TaskCompleted)
        fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
    else
        fprintf(' %d', ck_Computed_bits);
        fprintf('\n');  
        error = sum(abs(max(Task_Bit_Vec - ck_Computed_bits, 0)));
        fprintf('No Pass! sum:%d, size:%d error:%d\n',sum(Res_Check_TaskCompleted(:)), length(Res_Check_TaskCompleted), error);
        if error >= Num_UAV*1
            Final_Check = Final_Check + 1;
        end
    end

    %% [8] Association: Constraint of TAU
    disp('[8] Association: Constraint of TAU [8.1 8.2]');
    Res_Check_TAU1 = (Given_TAU_umn * (Matrix_Replicate_4_40') <= 1);
    [tmp_num_row, tmp_num_col] = size(Res_Check_TAU1);
    if sum(Res_Check_TAU1(:)) == tmp_num_row * tmp_num_col
        fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
    else
        fprintf('sum:%d, size:%d %d\n',sum(Res_Check_TAU1(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

    Res_Check_TAU2 = (Given_TAU_umn * (Matrix_Replicate_10_40') <= 1);
    [tmp_num_row, tmp_num_col] = size(Res_Check_TAU2);
    if sum(Res_Check_TAU2(:)) == tmp_num_row * tmp_num_col
        fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
    else
        fprintf('sum:%d, size:%d %d\n',sum(Res_Check_TAU2(:)), tmp_num_row, tmp_num_col);
        Final_Check = Final_Check + 1;
    end

%     %% [9] Local: Computing Constraint of L
%     disp('[9] Local: Computing Constraint of L');
%     ck_sum_L = ones(1,N) * Given_L_un;
%     Res_Check_L = (ck_sum_L <= 1);
%     [tmp_num_row, tmp_num_col] = size(Res_Check_L);
%     if sum(Res_Check_L(:)) == tmp_num_row * tmp_num_col
%         fprintf('Pass....%d %d\n',tmp_num_row, tmp_num_col);
%     else
%         fprintf('sum:%d, size:%d %d\n',sum(Res_Check_L(:)), tmp_num_row, tmp_num_col);
%         ck_sum_L
%     end
end

