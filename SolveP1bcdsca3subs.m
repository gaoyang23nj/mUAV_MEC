% traj(delay)-P(delay)-tauL(delay), not decrease then break
% 1. Optimize Trajectory (minimize Weighted Delay), 
%    SCA multiple r, inner loop, not decrease or infeasible
% 2. Optimize P (minimize Weighted Delay), SCA multiple r, inner loop
%    SCA multiple r, inner loop, not decrease or infeasible
% 3. Optimize th\TAU and L (minimize Weighted Delay);


clear;
clc;

diary('outputlog-AoTAlg.txt');
diary on;

tic;


%% Params
T = 100;
delta = 0.5;
N = 200;
% altidue 100m
H = 100;
% max velocity, 30m/s
v_max = 30;
% min distance between two UAVs, 50m
d_min = 100;
Num_UAV = 4;
Num_User = 10;
MAX_X = 1000;
MAX_Y = 1000;
%computation parameter
c_u = 1e3;
% local 10s, 20 time-slots
CPUFreq_User = 1e9;
% CPUFreq_User = 5e8;
CPUFreq_UAV = 10e9;
kappa_user = 1e-27;
kappa_uav = 1e-27;
%communication parameter (1M Hz, 0.1W)
rho = 1e-6;
Sigma2 = 1e-14;
Pu_max = 0.1;
Bandwidth = 1e6;

%constraint
E_user_max = 200;
E_uav_OE_max = 500;
E_uav_prop_max = 5000;

% % !!! the random seed
rng_seed = 1000;
rng(rng_seed);
Lk = 1e7;
Task_Bit_Vec = ones(1,Num_User)*Lk;
Loc_User_x = rand(1,Num_User)*MAX_X
Loc_User_y = rand(1,Num_User)*MAX_Y

%% hovering Params
% Utip 120m/s
prop_param_Utip = 120;

% profile drag coefficient
prop_param_delta = 0.012;
% rho, air density = 1.225 kg/m^3
prop_param_rho = 1.225;
% s, rotor solidity, 0.05
prop_param_s = 0.05;
% A, rotor disc area
prop_param_A = 0.503;
% blade angular velocity in radians/second
prop_param_Omega = 300;
% Rator radius in m
prop_param_R = 0.4;
prop_param_P0 = (prop_param_delta / 8) * prop_param_rho * prop_param_s * prop_param_A * power(prop_param_Omega, 3) * power(prop_param_R, 3);

% incremental correction factor to induced power
prop_param_k = 0.1;
% W, aircraft weight in Newton
prop_param_W = 20;
prop_param_Pi = (1 + prop_param_k) * power(prop_param_W, 3 / 2) / sqrt(2 * prop_param_rho * prop_param_A);

% d0, fuselage drag ratio
prop_param_d0 = 0.6;
% mean rotor induced velocity
prop_param_v0 = 4.03;

%% Useful Matrix
Matrix_Replicate_10_40 = zeros(Num_User, Num_User * Num_UAV);
for u=1:Num_User
    for m=1:Num_UAV
        Matrix_Replicate_10_40(u, u+(m-1)*Num_User) = 1;
    end
end
Matrix_Replicate_4_40 = zeros(Num_UAV, Num_User * Num_UAV);
for m=1:Num_UAV
    Matrix_Replicate_4_40(m, (m-1)*Num_User+1:m*Num_User) = 1;
end
Matrix_delete_user = zeros(Num_User * Num_UAV, Num_User * Num_UAV);
for m=1:Num_UAV
    Matrix_delete_user((m-1)*Num_User+1:m*Num_User, (m-1)*Num_User+1:m*Num_User) = 1;
end
Matrix_delete_user = Matrix_delete_user - eye(Num_User * Num_UAV);

%% Given Value
Given_TAU_umn = ones(N, Num_User * Num_UAV) / Num_User;
Given_L_un = ones(N, Num_User) * min(1.5 / N, CPUFreq_User*delta /(Task_Bit_Vec(1)*c_u));

Given_P_un = ones(N, Num_User).* Pu_max;
Given_F_umn = ones(N, Num_User * Num_UAV) * CPUFreq_UAV / Num_User;

Given_Q_mn_x = zeros(N,Num_UAV);
Given_Q_mn_y = zeros(N,Num_UAV);
Given_Qinit_mn_x = zeros(1,Num_UAV);
Given_Qinit_mn_y = zeros(1,Num_UAV);
init_center_loc = [MAX_X*0.25, MAX_Y*0.75; MAX_X*0.25, MAX_Y*0.25; MAX_X*0.75, MAX_Y*0.75; MAX_X*0.75, MAX_Y*0.25];
init_r = MAX_X * 0.25;
for m=1:Num_UAV
    Given_Qinit_mn_x(1, m) = init_center_loc(m,1) + cos(0)*init_r;
    Given_Qinit_mn_y(1, m) = init_center_loc(m,2) + sin(0)*init_r;
    for i=1:N
        tmp_theta = 2 * pi * (i/(N+1));
        Given_Q_mn_x(i, m) = init_center_loc(m,1) + cos(tmp_theta)*init_r;
        Given_Q_mn_y(i, m) = init_center_loc(m,2) + sin(tmp_theta)*init_r;
    end
end

%% Main Solve Programm...
MAX_Iteration = 50;
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

% Record Trajectory vs. Iteration
% Figure (Given_q_mn_x Given_q_mn_y)
figure(1);
hold on;
c = ['b';'k';'m';'r'];
for m=1:Num_UAV
    scatter(Given_Q_mn_x(:,m), Given_Q_mn_y(:,m),'.');
    % scatter(Given_Qinit_mn_x(:,m), Given_Qinit_mn_y(:,m),'*');
    scatter(Given_Qinit_mn_x(:,m), Given_Qinit_mn_y(:,m),'pentagram','LineWidth',3,'MarkerEdgeColor',c(m,:));
end
for u=1:Num_User
    % scatter(Loc_User_x(:,u), Loc_User_y(:,u),'^');
    scatter(Loc_User_x(:,u), Loc_User_y(:,u), '^','LineWidth',2);
end
xlabel('x(m)');
ylabel('y(m)');
box on;
set(gca, 'Fontname', 'Times New Roman','FontSize',12);


cvx_solver Mosek_2 

for iteration = 1:MAX_Iteration
    

    %% [2] SCA and Optimize Trajectory  ...
    Max_iteration_traj = 50;
    Record_target_traj = ones(1,Max_iteration_traj)*-inf;
    for iteration_traj=1:Max_iteration_traj
       %% [2.1] Q_r Y_r in SCA method for CVX Trajectory
        Given_Q_mn_x_r = Given_Q_mn_x;
        Given_Q_mn_y_r = Given_Q_mn_y;
        % init Y_slack_r
        dist_nplus1_x_r = [Given_Q_mn_x_r; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_x_r];
        dist_nplus1_y_r = [Given_Q_mn_y_r; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_y_r];
        Yslack_r_velocity = sqrt(power(dist_nplus1_x_r,2)+power(dist_nplus1_y_r,2))/delta;
        Yslack_r_term1 = sqrt(1 + 0.25 * power(Yslack_r_velocity/prop_param_v0, 4));
        Yslack_r_term2 = 0.5 * power(Yslack_r_velocity / prop_param_v0, 2);
        Given_Y_Slack_r = sqrt(Yslack_r_term1 - Yslack_r_term2);
        %Given_Y_Slack_r = zeros(N+1, Num_UAV);

       %% [2.2] CVX Trajectory (without bits constraint)
        disp('CVX Trajectory');
        cvx_begin
        %     cvx_solver SDPT3
            %cvx_solver SeDuMi
    %         cvx_solver_settings('write', 'dump.task.gz');
        %     cvx_solver_settings('MSK_IPAR_INFEAS_REPORT_AUTO', 'MSK_ON');
            cvx_precision low
        %     cvx_solver_settings('MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-2);

            variable Var_Q_mn_x(N, Num_UAV) nonnegative
            variable Var_Q_mn_y(N, Num_UAV) nonnegative
            variable Var_Y_Slack(N+1, Num_UAV) nonnegative
%             variable Var_S_Slack(N, Num_UAV*Num_User) nonnegative
%             variable Var_Ctl_Slack(N, Num_UAV*Num_User) nonnegative

            Matrix_sum_xy = zeros(2*Num_UAV, Num_UAV);
            for m=1:Num_UAV
                Matrix_sum_xy((m-1)*2+1:m*2, m) = 1;
            end

            % [1] UAV propulsion energy constraint
            % 1.1 slack with energy constraint
            dist_nplus1_x = [Var_Q_mn_x; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Var_Q_mn_x];
            dist_nplus1_y = [Var_Q_mn_y; Given_Qinit_mn_y] - [Given_Qinit_mn_y; Var_Q_mn_y];
            dist_power2 = pow_abs(dist_nplus1_x, 2) + pow_abs(dist_nplus1_y, 2);
            prop_energy_term1 = (1 + dist_power2*3)/(delta * delta * prop_param_Utip * prop_param_Utip) * prop_param_P0;
            dist_power3 = pow_pos(dist_power2, 1.5);
            prop_energy_term2 = dist_power3/(delta*delta*delta) * prop_param_d0 * prop_param_rho * prop_param_A * prop_param_s * 0.5;
            prop_energy_term3 = Var_Y_Slack * prop_param_Pi;
            E_uav_prop = sum(prop_energy_term1 + prop_energy_term2 + prop_energy_term3, 1) * delta;
            % E_uav_prop <= E_uav_prop_max;
            % [1.2] Y_slack constraint
            Y_term1 = power(Given_Y_Slack_r, 2) + 2 * Given_Y_Slack_r .* (Var_Y_Slack - Given_Y_Slack_r);
            dist_nplus1_x_r = [Given_Q_mn_x_r; Given_Qinit_mn_x] - [Given_Qinit_mn_x; Given_Q_mn_x_r];
            dist_nplus1_y_r = [Given_Q_mn_y_r; Given_Qinit_mn_y] - [Given_Qinit_mn_y; Given_Q_mn_y_r];
            dist_power2_r = power(dist_nplus1_x_r, 2) + power(dist_nplus1_y_r, 2); 
            Y_term2_1 = - dist_power2_r/(delta * delta * prop_param_v0 * prop_param_v0);
            Y_term2_2 = ((dist_nplus1_x_r .* dist_nplus1_x) + (dist_nplus1_y_r .* dist_nplus1_y)) * 2 / power(prop_param_delta*prop_param_v0, 2);
            %Con_Yslack = Y_term1 + Y_term2_1 + Y_term2_2 - pow_p(Var_Y_Slack, -2);
            Con_Yslack = Y_term1 + Y_term2_1 + Y_term2_2;
            %Con_Yslack >= pow_p(Var_Y_Slack, -2);

            % [2] offloading energy
            dist_uav_user_x_r = Given_Q_mn_x_r * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
            dist_uav_user_y_r = Given_Q_mn_y_r * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
            denomin_r = power(dist_uav_user_x_r, 2) + power(dist_uav_user_y_r, 2) + H*H;
            % Rate hat Taylor
            Rate_hat_Taylor_lg = (((Given_P_un * rho * Matrix_Replicate_10_40 / Sigma2) ./ denomin_r  ) * Matrix_Replicate_4_40') + 1;
            Rate_hat_Taylor_term1 = log2(Rate_hat_Taylor_lg) + log2(Sigma2);
            Rate_hat_Taylor_mul = ((Given_P_un * rho * Matrix_Replicate_10_40 * log2(exp(1)) / Sigma2) ./ power(denomin_r, 2)) ./ (Rate_hat_Taylor_lg * Matrix_Replicate_4_40);
            %dist_uav_user_x = Var_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
            %dist_uav_user_y = Var_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
            dist_uav_uav_x_r = Var_Q_mn_x * Matrix_Replicate_4_40 - Given_Q_mn_x_r * Matrix_Replicate_4_40;
            dist_uav_uav_y_r = Var_Q_mn_y * Matrix_Replicate_4_40 - Given_Q_mn_y_r * Matrix_Replicate_4_40;
            Rate_hat_Taylor_mul_dist = ((dist_uav_user_x_r .* dist_uav_uav_x_r) + (dist_uav_user_y_r .* dist_uav_uav_y_r))*2;
            Rate_hat_Taylor_term2 = (Rate_hat_Taylor_mul .* Rate_hat_Taylor_mul_dist) * Matrix_Replicate_4_40';
            Rate_hat_Taylor = Rate_hat_Taylor_term1 - Rate_hat_Taylor_term2;
            % Rate tilde Taylor
            Rate_tilde_Taylor_lg =  (((Given_P_un * rho * Matrix_Replicate_10_40 / Sigma2) ./ denomin_r) * Matrix_delete_user) + 1;
            Rate_tilde_Taylor_term1 = log2(Rate_tilde_Taylor_lg) + log2(Sigma2);
            Rate_tilde_Taylor_mul = ((Given_P_un * rho * Matrix_Replicate_10_40 * log2(exp(1)) / Sigma2) ./ power(denomin_r, 2)) ./ Rate_tilde_Taylor_lg;
            Rate_tilde_Taylor_term2 = (Rate_tilde_Taylor_mul .* Rate_hat_Taylor_mul_dist) * Matrix_delete_user;
            Rate_tilde_Taylor = Rate_tilde_Taylor_term1 - Rate_tilde_Taylor_term2;
            % Rate Taylor without Bandwidth
        %     Rate_comp_2taylor = Rate_hat_Taylor * Matrix_Replicate_4_40 - Rate_tilde_Taylor;
            Rate_comp_2taylor_true = Rate_hat_Taylor * Matrix_Replicate_4_40 - Rate_tilde_Taylor;
    %         e_comp_energy_ele = (((kappa_uav * c_u * delta * power(Given_F_umn, 2) * Bandwidth) .* Given_TAU_umn)  .* Rate_comp_2taylor) * Matrix_Replicate_4_40';
            % calculate the upper bound of Rate; if the upper bound of Rate can statisfy the user energy constraint, the real Rate also satisfies user energy constraint.
            e_comp_energy_ele = kappa_uav * power(Bandwidth, 3) * power(c_u,3) * delta * Given_TAU_umn .* pow_pos(Rate_comp_2taylor_true, 3);
            E_uav_comp = sum(e_comp_energy_ele, 1);
            % E_uav_comp <= E_uav_OE_max;

        %     tmp1 = (Given_P_un * rho / Sigma2 * Matrix_Replicate_10_40) .* inv_pos(Var_S_Slack + power(H,2));
        %     Rate_tilde_Sslack_lg = (tmp1 * Matrix_delete_user) + 1;
        %     Rate_tilde_Sslack = log2(Rate_tilde_Sslack_lg) + log2(Sigma2);
        %     Rate_comp_1taylor =  Rate_hat_Taylor * Matrix_Replicate_4_40 - Rate_tilde_Sslack;
        %     tmp1 = (Given_P_un * rho / Sigma2 * Matrix_Replicate_10_40) .* inv_pos(Var_S_Slack + power(H,2) );
        %     Rate_tilde_Sslack_lg = (tmp1 * Matrix_delete_user) + 1;
        %     Rate_tilde_Sslack = rel_entr(1, Rate_tilde_Sslack_lg)/log(2) - log2(Sigma2);
        % 
        %     tmp1 = (Given_P_un * rho / Sigma2 * Matrix_Replicate_10_40) .* inv_pos(Var_S_Slack + power(H,2) );
        %     Rate_tilde_Sslack_lg = inv_pos(tmp1 * Matrix_delete_user) + 1;
        %     Rate_tilde_Sslack = rel_entr(1, Rate_tilde_Sslack_lg)/log(2) - log2(Sigma2);
        %     Rate_comp_1taylor =  Rate_hat_Taylor * Matrix_Replicate_4_40 + Rate_tilde_Sslack;

            % [3] offloading capacity constraint
            Rate_comp_2taylor = max(Rate_comp_2taylor_true, 0);
            % Rate_comp_2taylor <= CPUFreq_UAV/Bd/c_u;

            % [4] mobility: maximum velocity constraint
            Dist_indelta = dist_power2;
            % Dist_indelta <= power(delta*v_max, 2)
            Vel_indelta = dist_power2 / power(delta*v_max, 2);
            % Vel_indelta <= 1;

            % [5] mobility: minimum distance between 2 UAVs constraint 
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
            St_2UAV_term1 = -(power(Given_Q_mn_x_r * Matrix_2UAV, 2) + power(Given_Q_mn_y_r * Matrix_2UAV, 2));
            St_2UAV_term2 = ((Given_Q_mn_x_r * Matrix_2UAV) .* (Var_Q_mn_x * Matrix_2UAV)) + ((Given_Q_mn_y_r * Matrix_2UAV) .* (Var_Q_mn_y * Matrix_2UAV));
        %     St_2UAV = St_2UAV_term1 + St_2UAV_term2 * 2;
            % St_2UAV >= d_min * d_min
            % Note!!! Scaling 1e5; [Rang in 1e-3 ~ 1e3] to solve numerically
            St_2UAV = (St_2UAV_term1/1e5) + (St_2UAV_term2 * 2/1e5);
            % St_2UAV >= (d_min * d_min)/1e5;

            % [6] task completed constraint 
            stcomp_offload_bits = Rate_comp_2taylor_true .* Given_TAU_umn * delta * (Matrix_Replicate_10_40') * Bandwidth /1e6;
            Computed_bits = sum(stcomp_offload_bits + Given_L_un * diag(Task_Bit_Vec/1e6), 1);
%             maximize(min(Computed_bits))

            % Optimization Target -- minmize the maximum weighted delay
%             target_offload_bits = Rate_comp_2taylor .* Given_TAU_umn * delta * (Matrix_Replicate_10_40') * Bandwidth /1e6;
%             target_time_bits = (target_offload_bits ./ repmat(Task_Bit_Vec/1e6,[N,1])) + Given_L_un;
%             Delay_Utility = sum(diag(1:1:N) * target_time_bits, 1);
%             Target = max(Delay_Utility);
%             minimize( Target );
            
            target_offload_bits = Rate_comp_2taylor .* Given_TAU_umn * delta * (Matrix_Replicate_10_40') * Bandwidth /1e6;
            target_time_bits = (target_offload_bits ./ repmat(Task_Bit_Vec/1e6,[N,1])) + Given_L_un;
            Delay_Utility = sum(diag(1:1:N) * target_time_bits, 1);
            Target = max(Delay_Utility);
            minimize( Target );
            
            subject to
                Var_Q_mn_x <= MAX_X
                Var_Q_mn_y <= MAX_Y
                E_uav_prop <= E_uav_prop_max
                Var_Y_Slack >= 0
                % important sub-area
                max(pow_abs(Var_Q_mn_x - Given_Q_mn_x_r, 2) + pow_abs(Var_Q_mn_y - Given_Q_mn_y_r, 2)) <= 100
%                 max(pow_abs(Var_Y_Slack - Given_Y_Slack_r, 2)) <= 0.1
                Con_Yslack >= pow_p(Var_Y_Slack, -2)
                E_uav_comp <= E_uav_OE_max
        %         Rate_comp_2taylor_true >= 0
                Rate_comp_2taylor <= CPUFreq_UAV/Bandwidth/c_u
                Vel_indelta <= 1
                St_2UAV >= (d_min * d_min)/1e5
                Computed_bits >= (Task_Bit_Vec/1e6)
        cvx_end
%         cvx_status
%         cvx_optval
       %% [2.3] Record after CVX tarjectory
        fprintf('OptRes_Traj [%d]-[%d]:%s %f\n', iteration, iteration_traj, cvx_status, cvx_optval);
        if isequal(cvx_status, 'Solved')
            Record_target_traj(1, iteration_traj) = cvx_optval;
            % tmp_change_target_traj = Record_target_traj(1, iteration_traj) - Record_target_traj(1, iteration_traj-1)
            % if converge then break (delay should not increase)
            %if ((iteration_traj >=2) && !( Record_target_traj(1, iteration_traj) <= Record_target_traj(1, iteration_traj-1)-0.001 ) )
            if ((iteration_traj >=2) && ( Record_target_traj(1, iteration_traj) > Record_target_traj(1, iteration_traj-1)-0.001 ) )
                fprintf('inner loop Traj Break! no decrement\n');
                break
            end
            
            % Record Trajectory vs. Iteration
            % Figure (Given_q_mn_x Given_q_mn_y)
            figure(2);
            hold on;
            c = ['b';'k';'m';'r'];
            for m=1:Num_UAV
                scatter(Var_Q_mn_x(:,m), Var_Q_mn_y(:,m),'.');
                % scatter(Given_Qinit_mn_x(:,m), Given_Qinit_mn_y(:,m),'*');
                scatter(Given_Qinit_mn_x(:,m), Given_Qinit_mn_y(:,m),'pentagram','LineWidth',3,'MarkerEdgeColor',c(m,:));
            end
            for u=1:Num_User
                % scatter(Loc_User_x(:,u), Loc_User_y(:,u),'^');
                scatter(Loc_User_x(:,u), Loc_User_y(:,u), '^','LineWidth',2);
            end
            xlabel('x(m)');
            ylabel('y(m)');
            box on;
            set(gca, 'Fontname', 'Times New Roman','FontSize',12);



            Given_Q_mn_x = Var_Q_mn_x;
            Given_Q_mn_y = Var_Q_mn_y;

            Record_allRes(2, iteration) = cvx_optval;
        else
            % infeasible
            fprintf('inner loop Traj Break! infeasible\n');
            ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
            [result,ck_Delay_Utility,ck_real_Delay_Utility,ck_prop_offload] = GetTargetValue(ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, N, Num_User, Num_UAV);
            CheckProc_func(Num_User, Num_UAV,ck_Rate*Bandwidth, Given_TAU_umn,Given_L_un,Given_P_un,Given_Q_mn_x,Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
            fprintf('ck_Rate_2taylor_true Analyzing... ...\n');
            ck_Rate_2taylor_true = Get2TaylorRate(Given_Q_mn_x, Given_Q_mn_y, Given_Q_mn_y, Given_Q_mn_x_r, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
            [result,ck_Delay_Utility,ck_real_Delay_Utility,ck_prop_offload] = GetTargetValue(ck_Rate_2taylor_true*Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, N, Num_User, Num_UAV);
            CheckProc_func(Num_User, Num_UAV,ck_Rate_2taylor_true*Bandwidth, Given_TAU_umn,Given_L_un,Given_P_un,Given_Q_mn_x,Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
            break
        end
        

    end

    
    %% [3] SCA and Optimize P...
    Max_iteration_P = 50;
    Record_target_P = ones(1,Max_iteration_traj)*-inf;
    for iteration_P=1:Max_iteration_P
        %% [3.1] P_r SCA method for CVX TAU and L
        Given_P_un_r = Given_P_un;

        %% [3.2] CVX P and induce F
        disp('CVX P and induce F');
        cvx_begin
            %cvx_precision low
            variable Var_P_un(N, Num_User) nonnegative
            cvx_precision low
        %     mosek_params = {'MSK_DPAR_INTPNT_CO_TOL_REL_GAP':  1e-2},
            cvx_solver_settings('MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-2);

            %[1] local computing energy
            le = power(Task_Bit_Vec * c_u,3) * kappa_user / (delta * delta);
            energy_user_loc_comp = power(Given_L_un, 3) * diag(le);
            energy_user_comm = Given_TAU_umn * (Matrix_Replicate_10_40') .* Var_P_un * delta;
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
            st_comp_energy_ele = (kappa_uav * power(Bandwidth,3) * delta  * power(c_u, 3) * Given_TAU_umn ) .* pow_pos(st_comm_Rate, 3);
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
            offload_bits = (st_Rate_concave .* Given_TAU_umn) * delta * (Matrix_Replicate_10_40') * (Bandwidth/1e6);
            local_bits = Given_L_un * diag(Task_Bit_Vec)/1e6;
            Computed_bits = sum(offload_bits + local_bits, 1);
            % Computed_bits >= Task_Bit_Vec
%             maximize( min(Computed_bits) )
            
        %     % Target
            target_comp_Rate = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde;
        %    target_comp_Rate = (Rate_hat_Taylor * Matrix_Replicate_4_40) - Rate_tilde_Taylor;
        %     % unit is Mbps
            target_offload_bits = (target_comp_Rate .* Given_TAU_umn) * delta * Matrix_Replicate_10_40' * Bandwidth/1e6;
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
            Given_P_un = Var_P_un;
            Record_allRes(3, iteration) = cvx_optval;
        else
            fprintf('inner loop P Break! infeasible\n');
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
    cvx_begin
        %cvx_precision low
        variable Var_L_un(N, Num_User) nonnegative
        variable Var_TAU_umn(N, Num_UAV*Num_User) nonnegative

        %[1] local computing energy
%         energy_user_loc_comp = Var_L_un * diag(power(Task_Bit_Vec * c_u, 3) * kappa_user /(delta*delta));
        energy_user_loc_comp = pow_pos(Var_L_un,3) * diag(power(Task_Bit_Vec * c_u, 3) * kappa_user /(delta*delta));
        energy_user_comm = Var_TAU_umn * (Matrix_Replicate_10_40') .* Given_P_un * delta;
        E_user = sum(energy_user_loc_comp + energy_user_comm,1);

        %[2] offloading computing energy
        energy_off_comp = power(Rate,3) * kappa_uav * power(c_u, 3) * delta .* Var_TAU_umn * (Matrix_Replicate_4_40');
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
        %Var_L_un * diag(Task_Bit_Vec * c_u) <= CPUFreq_User * delta

        %[5] task finished constraint
        offload_bits = Rate .* Var_TAU_umn * delta * (Matrix_Replicate_10_40')/1e6;
        local_bits = Var_L_un * diag(Task_Bit_Vec)/1e6;
        tmp = offload_bits + local_bits;
        Computed_bits = sum(offload_bits + local_bits, 1);
        %Computed_bits >= Task_Bit_Vec

        % Optimization Target
        offload_bits = Rate .* Var_TAU_umn * delta * (Matrix_Replicate_10_40')/1e6;
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
            Var_L_un * diag(Task_Bit_Vec/1e6 * c_u) <= CPUFreq_User/1e6 * delta
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
    end
    
    %% [1.4] Record Value after CVX TAU and L 
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
    
    % Complement Bits (modify Given_L_un)
    % Accurate Rate
    ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
    Given_L_un = ProcessL(ck_Rate, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, CPUFreq_User, c_u, N, Num_User, Num_UAV, Bandwidth);
    % Get Value
    [result,ck_Delay_Utility,ck_real_Delay_Utility,ck_prop_offload] = GetTargetValue(ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, N, Num_User, Num_UAV);
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

% Figure (Given_q_mn_x Given_q_mn_y)
% the final result
figure(3);
hold on;
title('final');
for m=1:Num_UAV
    scatter(Given_Q_mn_x(:,m), Given_Q_mn_y(:,m),'.');
    scatter(Given_Qinit_mn_x(:,m), Given_Qinit_mn_y(:,m),'*');
end
for u=1:Num_User
    scatter(Loc_User_x(:,u), Loc_User_y(:,u),'^');
end
box on;
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

% Figure (Given_q_mn_x Given_q_mn_y)
% the best result
figure(4);
hold on;
title('best');
for m=1:Num_UAV
    scatter(Record_min_Given_Q_mn_x(:,m), Record_min_Given_Q_mn_y(:,m),'.');
    scatter(Given_Qinit_mn_x(:,m), Given_Qinit_mn_y(:,m),'*');
end
for u=1:Num_User
    scatter(Loc_User_x(:,u), Loc_User_y(:,u),'^');
end
box on;
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

%Figure
figure(5);
hold on;
%plot(1:MAX_Iteration, Record_Res_iteration,'Marker','+','Color','b');
plot(1:MAX_Iteration, Record_Res_real_iteration,'Marker','o','Color','k');
xlabel('Iteration');
ylabel('Max-min AoT')
legend('cvxoptval','real value (after L)')
box on;
set(gca, 'Fontname', 'Times New Roman','FontSize',12);

save('111.mat');

Record_allRes;
CheckProc_func(Num_User, Num_UAV,ck_Rate*Bandwidth, Given_TAU_umn,Given_L_un,Given_P_un,Given_Q_mn_x,Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
            
toc;

diary off;
% filename = sprintf('alldatarngseed_rngseed%d.mat', rng_seed);
% save(filename)
