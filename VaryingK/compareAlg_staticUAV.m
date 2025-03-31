%% Static UAV H=100m

%% Main Solve Programm...
%     MAX_Iteration = 100;
Record_allRes = ones(3, MAX_Iteration) * (-1);
% record the output of Sub-Problem-1
Record_Res_iteration = ones(1,MAX_Iteration) * (-1);
Record_min_iteration = -1;
Record_min_cvxoptval = inf;
Record_min_Given_Q_mn_x = zeros(N,Num_UAV);
Record_min_Given_Q_mn_y = zeros(N,Num_UAV);
Record_min_TAU_umn = ones(N, Num_User * Num_UAV)* -1;
Record_min_L_un = ones(N, Num_User) * -1;
Record_min_P_un = ones(N, Num_User) * -1;

% record the 'completed-bits' output of Sub-Problem-1
Record_Res_real_iteration = ones(1,MAX_Iteration) * (-1);
Record_min_real_iteration = -1;
Record_min_result = inf;
Record_min_real_Given_Q_mn_x = zeros(N,Num_UAV);
Record_min_real_Given_Q_mn_y = zeros(N,Num_UAV);
Record_min_real_TAU_umn = ones(N, Num_User * Num_UAV)* -1;
Record_min_real_L_un = ones(N, Num_User) * -1;
Record_min_real_P_un = ones(N, Num_User) * -1;

for iteration = 1:MAX_Iteration

    %% [3] SCA and Optimize P...
    Max_iteration_P = 50;
    Record_target_P = ones(1,Max_iteration_traj)*-inf;
    for iteration_P=1:Max_iteration_P
        %% [3.1] P_r SCA method for CVX TAU and L
        Given_P_un_r = Given_P_un;

        %% [3.2] CVX P and induce F
        disp('CVX P and induce F');
        cvx_begin quiet
            %cvx_precision low
            variable Var_P_un(N, Num_User) nonnegative
            cvx_precision low
        %     mosek_params = {'MSK_DPAR_INTPNT_CO_TOL_REL_GAP':  1e-2},
            cvx_solver_settings('MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-2);

            %[1] local computing energy
            le = power(Task_Bit_Vec * c_u,3) * kappa_user / (Delta * Delta);
            energy_user_loc_comp = power(Given_L_un, 3) * diag(le);
            energy_user_comm = Given_TAU_umn * (Matrix_Replicate_10_40') .* Var_P_un * Delta;
            E_user = sum(energy_user_loc_comp + energy_user_comm,1);

            %[2] offloading computing energy
            % Rate hat
            dist_uav_user_x = Given_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
            dist_uav_user_y = Given_Q_mn_y * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
            denomin = power(dist_uav_user_x, 2) + power(dist_uav_user_y, 2) + H*H;
            Rate_hat_lg = ((Var_P_un * rho * Matrix_Replicate_10_40) /Sigma2 ./ denomin) * Matrix_Replicate_4_40' + 1;
            Rate_hat = log(Rate_hat_lg)/log(2) + log2(Sigma2);
        %     Rate_hat = -rel_entr(1,Rate_hat_lg)/log(2);
            % Rate tilde
            Rate_tilde_lg = ((Var_P_un * rho * Matrix_Replicate_10_40) /Sigma2 ./ denomin) * Matrix_delete_user + 1;
            Rate_tilde = log(Rate_tilde_lg) / log(2) + log2(Sigma2);
        %     Rate_tilde = -rel_entr(1,Rate_tilde_lg)/log(2);
            % Rate hat Taylor
            Rate_hat_Taylor_lg = ((Given_P_un_r * rho * Matrix_Replicate_10_40) / Sigma2 ./ denomin) * Matrix_Replicate_4_40' + 1;
            Rate_hat_Taylor_term1 = log2(Rate_hat_Taylor_lg) + log2(Sigma2);
            Rate_hat_Taylor_mul = (ones(N, Num_User) * Matrix_Replicate_10_40 * rho * log2(exp(1)) / Sigma2 ./ denomin) ./ (Rate_hat_Taylor_lg * Matrix_Replicate_4_40);
            Rate_hat_Taylor_term2 = Rate_hat_Taylor_mul .* ((Var_P_un - Given_P_un_r ) * Matrix_Replicate_10_40);
            Rate_hat_Taylor = Rate_hat_Taylor_term1 + (Rate_hat_Taylor_term2 * Matrix_Replicate_4_40');
            % Rate tilde Taylor
            Rate_tilde_Taylor_lg = ((Given_P_un_r * rho * Matrix_Replicate_10_40) /Sigma2 ./ denomin) * Matrix_delete_user + 1;
            Rate_tilde_Taylor_term1 = log2(Rate_tilde_Taylor_lg) + log2(Sigma2);
            Rate_tilde_Taylor_mul = (ones(N, Num_User) * Matrix_Replicate_10_40 * rho * log2(exp(1)) / Sigma2 ./ denomin) ./ Rate_tilde_Taylor_lg;
            Rate_tilde_Taylor_term2 = Rate_tilde_Taylor_mul .* ((Var_P_un - Given_P_un_r ) * Matrix_Replicate_10_40);
            Rate_tilde_Taylor = Rate_tilde_Taylor_term1 + (Rate_tilde_Taylor_term2 * Matrix_delete_user);
            % convex Rate
            st_comm_Rate = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde;
            st_comp_energy_ele = (kappa_uav * power(Bandwidth,3) * Delta  * power(c_u, 3) * Given_TAU_umn ) .* pow_pos(st_comm_Rate, 3);
            E_uav_comp = sum(st_comp_energy_ele * Matrix_Replicate_4_40', 1);
            % E_uav_comp <= E_uav_OE_max;

            %[3] offloading constraint
            st_Rate_convex = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde;
            Freq_Needed = max(st_Rate_convex/1e6, 0) * c_u;
        %     Freq_Needed <= CPUFreq_UAV

            %[4] task finished constraint
            % convex Rate
            st_Rate_concave = (Rate_hat * Matrix_Replicate_4_40) - Rate_tilde_Taylor;
        %    st_Rate_concave = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde_Taylor;
        %     st_Rate_concave = (Rate_hat_Taylor_term1 * Matrix_Replicate_4_40) - Rate_tilde_Taylor_term1;
            % unit is Mb
            offload_bits = (st_Rate_concave .* Given_TAU_umn) * Delta * (Matrix_Replicate_10_40') * (Bandwidth/1e6);
            local_bits = Given_L_un * diag(Task_Bit_Vec)/1e6;
            Computed_bits = sum(offload_bits + local_bits, 1);
            % Computed_bits >= Task_Bit_Vec
%             maximize( min(Computed_bits) )
            
        %     % Target
            target_comp_Rate = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde;
        %    target_comp_Rate = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde_Taylor;
        %     % unit is Mbps
            target_offload_bits = (target_comp_Rate .* Given_TAU_umn) * Delta * Matrix_Replicate_10_40' * Bandwidth/1e6;
            targe_time_bits = (target_offload_bits./repmat(Task_Bit_Vec/1e6, [N, 1])) + Given_L_un;
            Delay_Utility = sum(diag(1:1:N) * targe_time_bits, 1);
            Target = max(Delay_Utility);
            minimize( Target )

            % sub-area constraint
            subarea = max(pow_abs(Var_P_un - Given_P_un_r, 2));

            subject to
                Var_P_un <= Pu_max
                subarea <= 0.01 * Pu_max
                E_user <= E_user_max
                E_uav_comp <= E_uav_OE_max
                Freq_Needed <= CPUFreq_UAV/1e6
                Computed_bits >= Task_Bit_Vec/1e6
        cvx_end
%         cvx_status
%         cvx_optval
       %% [3.3] Record after CVX P
        fprintf('OptRes_P [%d]-[%d]:%s %f\n', iteration, iteration_P, cvx_status, cvx_optval);
        if isequal(cvx_status, 'Solved')
            Record_target_P(1, iteration_P) = cvx_optval;
            fprintf('transmit power [%f]~[%f]\n', min(min(Var_P_un)), max(max(Var_P_un)));
            % tmp_change_target_P = Record_target_P(1, iteration_P) - Record_target_P(1, iteration_P-1)
            % if converge (less than not OK!)then break  (delay should not increase)
            %if ((iteration_P >=2) && !( Record_target_P(1, iteration_P) <= Record_target_P(1, iteration_P-1)-0.001 ) )
            if ((iteration_P >=2) && ( Record_target_P(1, iteration_P) > Record_target_P(1, iteration_P-1)-0.001 ) )
                fprintf('inner loop P Break! no decrement\n');
                break
            end

            ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Var_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
            tmpCVXP_Given_L_un = ProcessL(ck_Rate, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, CPUFreq_User, c_u, N, Num_User, Num_UAV, Bandwidth);
            [Final_Check, ck_prop_energy] = CheckProc_func(Num_User, Num_UAV,ck_Rate*Bandwidth, Given_TAU_umn, tmpCVXP_Given_L_un, Given_P_un,Given_Q_mn_x, Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
            if Final_Check == 0
                fprintf('Pass! CVX Power Check \n');
                % record the result
                Given_P_un = Var_P_un;
                Record_allRes(3, iteration) = cvx_optval;
            else
                fprintf('Break! NoPass CVX Power Check (%d)\n', Final_Check);
                break;
            end
        else
            fprintf('inner loop Power Break! infeasible\n');
            break
        end
        
    end
    
    %% [1] Optimize TAU and L ...
    %% [1.1] accurate Rate for  CVX TAU and L
    % Rate Hat
    uav_x = Given_Q_mn_x * Matrix_Replicate_4_40;
    uav_y = Given_Q_mn_y * Matrix_Replicate_4_40;
    user_x = ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
    user_y = ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
    denomin = power(uav_x - user_x, 2) + power(uav_y - user_y, 2) + H*H;
    Rate_hat_lg = ((Given_P_un * rho * Matrix_Replicate_10_40) ./ denomin) * Matrix_Replicate_4_40' + Sigma2;
    Rate_hat = log2(Rate_hat_lg) .* Bandwidth;
    % Rate tilde
    Rate_tilde_lg = ((Given_P_un .* rho * Matrix_Replicate_10_40) ./ denomin) * Matrix_delete_user + Sigma2;
    Rate_tilde = log2(Rate_tilde_lg) .* Bandwidth;
    Rate = Rate_hat * Matrix_Replicate_4_40 - Rate_tilde;

    %% [1.2] CVX TAU and L
    disp('CVX TAU and L');
    cvx_begin quiet
        %cvx_precision low
        variable Var_L_un(N, Num_User) nonnegative
        variable Var_TAU_umn(N, Num_UAV*Num_User) nonnegative

        %[1] local computing energy
%         energy_user_loc_comp = Var_L_un * diag(power(Task_Bit_Vec * c_u, 3) * kappa_user /(Delta*Delta));
        energy_user_loc_comp = pow_pos(Var_L_un,3) * diag(power(Task_Bit_Vec * c_u, 3) * kappa_user /(Delta*Delta));
        energy_user_comm = Var_TAU_umn * (Matrix_Replicate_10_40') .* Given_P_un * Delta;
        E_user = sum(energy_user_loc_comp + energy_user_comm,1);

        %[2] offloading computing energy
        energy_off_comp = power(Rate,3) * kappa_uav * power(c_u, 3) * Delta .* Var_TAU_umn * (Matrix_Replicate_4_40');
        E_uav_comp = sum(energy_off_comp, 1);
        %[3] association constraint
        %[3.1]
        % 0 <= Var_TAU_umn <= 1
        %[3.2]
        %Var_TAU_umn * (Matrix_Replicate_4_40') <= 1
        %[3.3]
        %Var_TAU_umn * (Matrix_Replicate_10_40') <= 1

        %[4] offloading constraint
        %[4.1]
        %0 <= Var_L_un <= 1
        %[4.2]
        % ones(1,N) * Var_L_un <= 1
        %[4.3] cannot beyond the compuataion capacity
        %Var_L_un * diag(Task_Bit_Vec * c_u) <= CPUFreq_User * Delta

        %[5] task finished constraint
        offload_bits = Rate .* Var_TAU_umn * Delta * (Matrix_Replicate_10_40')/1e6;
        local_bits = Var_L_un * diag(Task_Bit_Vec)/1e6;
        tmp = offload_bits + local_bits;
        Computed_bits = sum(offload_bits + local_bits, 1);
        %Computed_bits >= Task_Bit_Vec

        % Optimization Target
        offload_bits = Rate .* Var_TAU_umn * Delta * (Matrix_Replicate_10_40')/1e6;
        offload_ratio = offload_bits ./ repmat(Task_Bit_Vec/1e6,[N,1]);
        Delay_Utility = sum(diag(1:1:N) * (offload_ratio + Var_L_un),1);
        Target = max(Delay_Utility);

        minimize( Target )
        subject to
            0 <= Var_L_un <= 1
            0 <= Var_TAU_umn <= 1
            E_user <= E_user_max
            E_uav_comp <= E_uav_OE_max
            Var_TAU_umn * (Matrix_Replicate_4_40') <= 1
            Var_TAU_umn * (Matrix_Replicate_10_40') <= 1
            ones(1,N) * Var_L_un <= 1
            Var_L_un * diag(Task_Bit_Vec/1e6 * c_u) <= CPUFreq_User/1e6 * Delta
            Computed_bits >= Task_Bit_Vec/1e6
    cvx_end
    fprintf('OptRes_TauL [%d]:%s %f\n', iteration, cvx_status, cvx_optval);
%     cvx_status
%     cvx_optval
    if isequal(cvx_status, 'Solved')
        Given_L_un = Var_L_un;
        Given_TAU_umn = Var_TAU_umn;
        num_trasf_tau = 0;
        num_trasf_L = 0;
        % let the very small TAU to be 0
        for i=1:N
            for j=1:(Num_UAV * Num_User)
                if Given_TAU_umn(i,j) <= 1e-5
                    Given_TAU_umn(i,j) = 0;
                    num_trasf_tau = num_trasf_tau + 1;
                end
            end
        end
        % let the very small L=0
        for i=1:N
            for j=1:Num_User
                if Given_L_un(i,j) <= 0
                    Given_L_un(i,j) = 0;
                    num_trasf_L = num_trasf_L + 1;
                end
            end
        end
        % Need to Arrange & Check the result 
        Record_allRes(1, iteration) = cvx_optval;
    else
        fprintf('ERROR! CVX Status!!! CVX TAU and L\n');
        break;
    end

    
    %% [1.4] Record Value after CVX TAU and L
    fprintf('CVX BCD-[%d]\n', iteration);
    % Complement Bits (modify Given_L_un)
    % Accurate Rate
    ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
    Given_L_un = ProcessL(ck_Rate, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, CPUFreq_User, c_u, N, Num_User, Num_UAV, Bandwidth);
    [Final_Check, ck_prop_energy] = CheckProc_func(Num_User, Num_UAV,ck_Rate*Bandwidth, Given_TAU_umn,Given_L_un,Given_P_un,Given_Q_mn_x,Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
    % Get Value
    [result,ck_Delay_Utility,ck_real_Delay_Utility,ck_prop_offload] = GetTargetValue(ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, N, Num_User, Num_UAV);
    
    if Final_Check == 0
        fprintf('Pass! Record Value after CVX TAU and L \n');
    else
        fprintf('Break! NoPass Record Value after CVX TAU and L (%d) \n', Final_Check);
        break;
    end

    Record_Res_iteration(1,iteration) = cvx_optval;
    % record the best one according to the Output of SubProblem1
    if cvx_optval < Record_min_cvxoptval
        Record_min_iteration = iteration;
        Record_min_cvxoptval = cvx_optval;
        Record_min_Given_Q_mn_x = Given_Q_mn_x;
        Record_min_Given_Q_mn_y = Given_Q_mn_y;
        Record_min_TAU_umn = Given_TAU_umn;
        Record_min_L_un = Given_L_un;
        Record_min_P_un = Given_P_un;
    end
    
    Record_Res_real_iteration(1,iteration) = result;
    if result < Record_min_result
        Record_min_real_iteration = iteration;
        Record_min_result = result;
        Record_min_real_Given_Q_mn_x = Given_Q_mn_x;
        Record_min_real_Given_Q_mn_y = Given_Q_mn_y;
        Record_min_real_TAU_umn = Given_TAU_umn;
        Record_min_real_L_un = Given_L_un;
        Record_min_real_P_un = Given_P_un;
    end

end

%Figure
figure(1);
hold on;
%plot(1:MAX_Iteration, Record_Res_iteration,'Marker','+','Color','b');
plot(1:MAX_Iteration, Record_Res_real_iteration,'Marker','o','Color','k');
xlabel('Iteration');
ylabel('Max-min AoT');
% legend('cvxoptval','real value (after L)');
legend('cvxoptval');
Record_allRes

A_finalIteration_Res = Record_Res_real_iteration(end);
A_breakiteration = iteration;
