% simulate
clear all;
clc;

diary('outputlog.txt');
diary on;

% choose Mosek aom
cvx_solver Mosek_2;

time_begin = clock;
fprintf('Now:')
fprintf(' %d',clock);
fprintf('\n');
tic;
Times = 2;

% rng seed
t = clock;
rng(t(6)*1000+t(5)*60+t(4)*3600);
fprintf('rng seed:%d,%f\n',(t(6)+t(5)*60+t(4)*3600)*1000, rand(1));

% conducting program
MAX_Iteration = 50;
Num_UAV = 4;
Num_init_K = 7;
Num_end_K = 14;
list_num_user = Num_init_K : Num_end_K;
L_k = 5e7;



%% record Result of varying K multiple times
% OPT Alg
all_Record_VaryingK_OPT =  ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_OPT_real = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_OPT_Math = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_OPT_Sim = ones(Times, length(list_num_user)) * (-1);
% Local Alg
all_Record_VaryingK_Local_Sim = ones(Times, length(list_num_user)) * (-1);
% Static Alg
all_Record_VaryingK_Static =  ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_Static_real = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_Static_Math = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_Static_Sim = ones(Times, length(list_num_user)) * (-1);
% Peak Power Alg
all_Record_VaryingK_PeakPower =  ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_PeakPower_real = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_PeakPower_Math = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_PeakPower_Sim = ones(Times, length(list_num_user)) * (-1);
% Upper Bound Cal.
all_Record_VaryingK_UpperBound = ones(Times, length(list_num_user)) * (-1);
% Lower Bound Cal.
all_Record_VaryingK_LowerBound = ones(Times, length(list_num_user)) * (-1);
% Maxmin OffBits Alg
all_Record_VaryingK_MaxminOffBits =  ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MaxminOffBits_real = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MaxminOffBits_Math = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MaxminOffBits_Sim = ones(Times, length(list_num_user)) * (-1);
% Min sum OPTime
all_Record_VaryingK_MinOPTime =  ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MinOPTime_real = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MinOPTime_Math = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MinOPTime_Sim = ones(Times, length(list_num_user)) * (-1);

% Total energy comsumed
all_Record_VaryingK_OPT_Energy = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_Local_Energy = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_Static_Energy = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_PeakPower_Energy = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MaxminOffBits_Energy = ones(Times, length(list_num_user)) * (-1);
all_Record_VaryingK_MinOPTime_Energy = ones(Times, length(list_num_user)) * (-1);


label_break = false;
for tmp_User = 1:length(list_num_user)
    if label_break
        break
    end
    fprintf('Num_User = %d\n', list_num_user(tmp_User));
    Num_User = list_num_user(tmp_User);
    Task_Bit_Vec = ones(1,Num_User)*L_k;
    INIT_PARAMS_K;
    for current_times = 1:Times

        Loc_User_x = rand(1,Num_User)*MAX_X;
        fprintf('Loc_User_x:');
        fprintf(' %d', Loc_User_x);
        Loc_User_y = rand(1,Num_User)*MAX_Y;
        fprintf('\nLoc_User_y:');
        fprintf(' %d', Loc_User_y);
        fprintf('\n');

        %% Local Alg
        disp('Local Alg');
        % how many bits in one slot
        bits_oneslot = CPUFreq_User * Delta / c_u;
        % Task_Bit_Vec / bits_oneslot
        remain_bits = Task_Bit_Vec;
        tmp_res_local = zeros(1, Num_User);
        tmp_energy = zeros(1, Num_User);
        for tmp_localcomp_user = 1:Num_User
            for i=1:N
                cal_bits = min(remain_bits(1, tmp_localcomp_user), bits_oneslot );
                tmp_res_local(1, tmp_localcomp_user) = tmp_res_local(1, tmp_localcomp_user) + cal_bits * i;
                remain_bits(1, tmp_localcomp_user) = remain_bits(1, tmp_localcomp_user) - cal_bits;
 				tmp_energy(1, tmp_localcomp_user) = tmp_energy(1, tmp_localcomp_user) + power(cal_bits * c_u, 3)*kappa_user/(Delta*Delta);
            end
        end
        tmp_res_local;
        tmp_res_local./Task_Bit_Vec;
        ResOPT = max(tmp_res_local);
        all_Record_VaryingK_Local_Sim(current_times, tmp_User) = ResOPT;
        Res_Total_Energy = sum(tmp_energy);
        all_Record_VaryingK_Local_Energy(current_times, tmp_User) = Res_Total_Energy;
        fprintf('LocalAlg ResOPT: %f %f(10^7); Total_Energy: %f\n', ResOPT, ResOPT/1e7, Res_Total_Energy);


        %% Lower Upper Bound of AoT
        % Upper Bound is Local Computing, (maybe with higher or worst Comm.)
        capacity = (CPUFreq_User * Delta)/c_u;
        fprintf('UP capac: %f %f(Mb)\n', capacity, capacity/1e6);
        remain_bits = Task_Bit_Vec;
        UB_AoT_k = zeros(1, Num_User);
        for idx_upB_usr = 1:Num_User
            tmpUpBound_N_up = ceil(Task_Bit_Vec(idx_upB_usr) / capacity);
            fprintf('tmpLwBound_N_lb: %f\n', tmpUpBound_N_up);
            for idx_slot = 1:N
                cal_bits = min(capacity, remain_bits(idx_upB_usr));
                UB_AoT_k(idx_upB_usr) = UB_AoT_k(idx_upB_usr) + idx_slot * cal_bits;
                remain_bits(idx_upB_usr) =  remain_bits(idx_upB_usr) - cal_bits;
            end
        end
        max(UB_AoT_k);
        fprintf('UpperBoundAoT ResLocal: %f %f(10^7)\n', max(UB_AoT_k), max(UB_AoT_k)/1e7);
        all_Record_VaryingK_UpperBound(current_times, tmp_User) = max(UB_AoT_k);    

        % Lower Bound with min(the best Comm., best Comp.)
        up_transmit_rate = Bandwidth * log2((Pu_max * rho) / (H*H)/Sigma2);
        up_process_rate = CPUFreq_UAV * Delta / c_u;
        UAV_capa = max(0, min(up_transmit_rate, up_process_rate))*Num_UAV/Num_User;
        % remote capacity + local capacity
        capacity = UAV_capa + (CPUFreq_User * Delta)/c_u;
        fprintf('LW capac: %f %f(Mb)\n', capacity, capacity/1e6);
        remain_bits = Task_Bit_Vec;
        LB_AoT_k = zeros(1, Num_User);
        for idx_lwB_usr = 1:Num_User
            tmpLwBound_N_lb = floor(Task_Bit_Vec(idx_lwB_usr) / capacity);
            fprintf('tmpLwBound_N_lb: %f\n', tmpLwBound_N_lb);
            for idx_slot = 1:N
                cal_bits = min(capacity, remain_bits(idx_lwB_usr));
                LB_AoT_k(idx_lwB_usr) = LB_AoT_k(idx_lwB_usr) + idx_slot * cal_bits;
                remain_bits(idx_lwB_usr) =  remain_bits(idx_lwB_usr) - cal_bits;
            end
        end
        min(LB_AoT_k);
        % list_lower_bound(index_num_user) = min(LB_AoT_k);
        fprintf('LowerBoundAoT ResLocal: %f %f(10^7)\n', min(LB_AoT_k), min(LB_AoT_k)/1e7);
        all_Record_VaryingK_LowerBound(current_times, tmp_User) = min(LB_AoT_k);


        %% OPT proposed Alg
        disp('OPT proposed Alg');
        % Given Value
        Given_TAU_umn = ones(N, Num_User * Num_UAV) / Num_User;
        Given_L_un = ones(N, Num_User) * min(1.5 / N, CPUFreq_User*Delta /(Task_Bit_Vec(1)*c_u));
        Given_P_un = ones(N, Num_User).* Pu_max;
        Given_F_umn = ones(N, Num_User * Num_UAV) * CPUFreq_UAV / Num_User;
        Given_Q_mn_x = zeros(N,Num_UAV);
        Given_Q_mn_y = zeros(N,Num_UAV);
        Given_Qinit_mn_x = zeros(1,Num_UAV);
        Given_Qinit_mn_y = zeros(1,Num_UAV);
        init_center_loc = [MAX_X*0.25, MAX_Y*0.75; MAX_X*0.25, MAX_Y*0.25; MAX_X*0.75, MAX_Y*0.75; MAX_X*0.75, MAX_Y*0.25];
        init_r = 250;
        for m=1:Num_UAV
            Given_Qinit_mn_x(1, m) = init_center_loc(m,1) + cos(0)*init_r;
            Given_Qinit_mn_y(1, m) = init_center_loc(m,2) + sin(0)*init_r;
            for i=1:N
                tmp_theta = 2 * pi * (i/(N+1));
                Given_Q_mn_x(i, m) = init_center_loc(m,1) + cos(tmp_theta)*init_r;
                Given_Q_mn_y(i, m) = init_center_loc(m,2) + sin(tmp_theta)*init_r;
            end
        end
        % conduct alg
        compareAlg_OPT;
        if A_breakiteration ~= MAX_Iteration
            fprintf('Break! compareAlg_OPT [K:%d]-[%d]', Num_User, current_times);
            label_break = true;
            break;
        end
        % Record Res
        all_Record_VaryingK_OPT(current_times, tmp_User) =  Record_Res_iteration(end);
        all_Record_VaryingK_OPT_real(current_times, tmp_User) = Record_Res_real_iteration(end);
        sprintf('OPTALG CAL_RES: %f %f', Record_Res_iteration(end), Record_Res_real_iteration(end));
        % Math Calculation
        ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
        [TargetStatic,Delay_UtilityStatic,real_Delay_UtilityStatic,prop_offloadStatic] = GetTargetValue(ck_Rate * Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, N, Num_User, Num_UAV);
        % currate Time with higher N
        fprintf('Static Math Calculation:%f', TargetStatic);
        all_Record_VaryingK_OPT_Math(current_times, tmp_User) = TargetStatic;
        % Simulate 
        assert(N == T/Delta);
        Task_Bit_Vec;
        offload_bits = ck_Rate * Bandwidth .* Given_TAU_umn * Delta * (Matrix_Replicate_10_40');
        % diag(1:1:N) * Delta * offload_bits
        local_bits = Given_L_un * diag(Task_Bit_Vec);
        % [local_Freq]
        total_bits = offload_bits + local_bits;
        sum(total_bits) >= Task_Bit_Vec;
        tmp_res_opt = zeros(1, Num_User);
        for i = 1:N
            tmp_res_opt = tmp_res_opt + total_bits(i,:) * i;
%             tmp_res_opt = tmp_res_opt + total_bits(i,:) * i * Delta;
        end
        tmp_res_opt
        tmp_res_opt./Task_Bit_Vec;
        ResOPT = max(tmp_res_opt);
        all_Record_VaryingK_OPT_Sim(current_times, tmp_User) = ResOPT;
        Res_Total_Energy = GetTotalEnergy(Num_User, Num_UAV, ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Given_P_un, Given_Q_mn_x, Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
        all_Record_VaryingK_OPT_Energy(current_times, tmp_User) = Res_Total_Energy;
        fprintf('OptAlg ResOPT: %f %f(10^7); Total_Energy: %f\n', ResOPT, ResOPT/1e7, Res_Total_Energy);


        %% Static UAV Alg
        disp('Static UAV Alg');
        % Given Value
        Given_TAU_umn = ones(N, Num_User * Num_UAV) / Num_User;
        Given_L_un = ones(N, Num_User) * min(1.5 / N, CPUFreq_User*Delta /(Task_Bit_Vec(1)*c_u));
        Given_P_un = ones(N, Num_User).* Pu_max;
        Given_F_umn = ones(N, Num_User * Num_UAV) * CPUFreq_UAV / Num_User;
        Given_Q_mn_x = zeros(N,Num_UAV);
        Given_Q_mn_y = zeros(N,Num_UAV);
        Given_Qinit_mn_x = zeros(1,Num_UAV);
        Given_Qinit_mn_y = zeros(1,Num_UAV);
        init_center_loc = [MAX_X*0.25, MAX_Y*0.75; MAX_X*0.25, MAX_Y*0.25; MAX_X*0.75, MAX_Y*0.75; MAX_X*0.75, MAX_Y*0.25];
        init_r = 250;
        for m=1:Num_UAV
            Given_Qinit_mn_x(1, m) = init_center_loc(m,1) + cos(0)*init_r;
            Given_Qinit_mn_y(1, m) = init_center_loc(m,2) + sin(0)*init_r;
            for i=1:N
                tmp_theta = 2 * pi * (i/(N+1));
                Given_Q_mn_x(i, m) = init_center_loc(m,1) + cos(tmp_theta)*init_r;
                Given_Q_mn_y(i, m) = init_center_loc(m,2) + sin(tmp_theta)*init_r;
            end
        end
        Given_Q_mn_x = repmat(Given_Qinit_mn_x, N,1);
        Given_Q_mn_y = repmat(Given_Qinit_mn_y, N,1);
        % conduct alg
        compareAlg_staticUAV;
        if A_breakiteration ~= MAX_Iteration
            fprintf('Break! compareAlg_staticUAV [K:%d]-[%d]', Num_User, current_times);
            label_break = true;
            break;
        end
        all_Record_VaryingK_Static(current_times, tmp_User) =  Record_Res_iteration(end);
        all_Record_VaryingK_Static_real(current_times, tmp_User) = Record_Res_real_iteration(end);
        sprintf('StaticUAVAlg CAL_RES: %f %f', Record_Res_iteration(end), Record_Res_real_iteration(end));
        % Math Calculation
        ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
        [TargetStatic,Delay_UtilityStatic,real_Delay_UtilityStatic,prop_offloadStatic] = GetTargetValue(ck_Rate * Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, N, Num_User, Num_UAV);
        % currate Time with higher N
        fprintf('StaticUAVAlg Math Calculation:%f', TargetStatic);
        all_Record_VaryingK_Static_Math(current_times, tmp_User) = TargetStatic;
        % Simulate 
        assert(N == T/Delta)
        Task_Bit_Vec
        offload_bitsStatic = ck_Rate * Bandwidth .* Given_TAU_umn * Delta * (Matrix_Replicate_10_40');
        % diag(1:1:N) * Delta * offload_bits
        local_bitsStatic = Given_L_un * diag(Task_Bit_Vec);
        % [local_Freq]
        total_bitsStatic = offload_bitsStatic + local_bitsStatic;
        sum(total_bitsStatic) >= Task_Bit_Vec;
        tmp_res_opt = zeros(1, Num_User);
        for i = 1:N
            tmp_res_opt = tmp_res_opt + total_bitsStatic(i,:) * i;
        %     tmp_res_opt = tmp_res_opt + total_bits(i,:) * i * Delta;
        end
        tmp_res_opt
        tmp_res_opt./Task_Bit_Vec;
        ResOPT = max(tmp_res_opt);
        all_Record_VaryingK_Static_Sim(current_times, tmp_User) = ResOPT;
        Res_Total_Energy = GetTotalEnergy(Num_User, Num_UAV, ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Given_P_un, Given_Q_mn_x, Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
        all_Record_VaryingK_Static_Energy(current_times, tmp_User) = Res_Total_Energy;
        fprintf('StaticUAV Alg Res: %f %f(10^7); Total_Energy: %f\n', ResOPT, ResOPT/1e7, Res_Total_Energy);


        %% Peak Power Alg
        disp('Peak Power Alg');
        % Given Value
        Given_TAU_umn = ones(N, Num_User * Num_UAV) / Num_User;
        Given_L_un = ones(N, Num_User) * min(1.5 / N, CPUFreq_User*Delta /(Task_Bit_Vec(1)*c_u));
        Given_P_un = ones(N, Num_User).* Pu_max;
        Given_F_umn = ones(N, Num_User * Num_UAV) * CPUFreq_UAV / Num_User;
        Given_Q_mn_x = zeros(N,Num_UAV);
        Given_Q_mn_y = zeros(N,Num_UAV);
        Given_Qinit_mn_x = zeros(1,Num_UAV);
        Given_Qinit_mn_y = zeros(1,Num_UAV);
        init_center_loc = [MAX_X*0.25, MAX_Y*0.75; MAX_X*0.25, MAX_Y*0.25; MAX_X*0.75, MAX_Y*0.75; MAX_X*0.75, MAX_Y*0.25];
        init_r = 250;
        for m=1:Num_UAV
            Given_Qinit_mn_x(1, m) = init_center_loc(m,1) + cos(0)*init_r;
            Given_Qinit_mn_y(1, m) = init_center_loc(m,2) + sin(0)*init_r;
            for i=1:N
                tmp_theta = 2 * pi * (i/(N+1));
                Given_Q_mn_x(i, m) = init_center_loc(m,1) + cos(tmp_theta)*init_r;
                Given_Q_mn_y(i, m) = init_center_loc(m,2) + sin(tmp_theta)*init_r;
            end
        end
        % conduct alg
        compareAlg_peakTranPower;
        if A_breakiteration ~= MAX_Iteration
            fprintf('Break! compareAlg_peakTranPower [K:%d]-[%d]', Num_User, current_times);
            label_break = true;
            break;
        end
        all_Record_VaryingK_PeakPower(current_times, tmp_User) =  Record_Res_iteration(end);
        all_Record_VaryingK_PeakPower_real(current_times, tmp_User) = Record_Res_real_iteration(end);
        sprintf('PeakPowerAlg CAL_RES: %f %f', Record_Res_iteration(end), Record_Res_real_iteration(end));
        % Math Calculation
        ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
        [TargetStatic,Delay_UtilityStatic,real_Delay_UtilityStatic,prop_offloadStatic] = GetTargetValue(ck_Rate * Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, N, Num_User, Num_UAV);
        % currate Time with higher N
        fprintf('PeakPowerAlg Math Calculation:%f', TargetStatic);
        all_Record_VaryingK_PeakPower_Math(current_times, tmp_User) = TargetStatic;
        % Simulate 
        assert(N == T/Delta)
        Task_Bit_Vec
        offload_bitsStatic = ck_Rate * Bandwidth .* Given_TAU_umn * Delta * (Matrix_Replicate_10_40');
        % diag(1:1:N) * Delta * offload_bits
        local_bitsStatic = Given_L_un * diag(Task_Bit_Vec);
        % [local_Freq]
        total_bitsStatic = offload_bitsStatic + local_bitsStatic;
        sum(total_bitsStatic) >= Task_Bit_Vec;
        tmp_res_opt = zeros(1, Num_User);
        for i = 1:N
            tmp_res_opt = tmp_res_opt + total_bitsStatic(i,:) * i;
        %     tmp_res_opt = tmp_res_opt + total_bits(i,:) * i * Delta;
        end
        tmp_res_opt
        tmp_res_opt./Task_Bit_Vec;
        ResOPT = max(tmp_res_opt);
        all_Record_VaryingK_PeakPower_Sim(current_times, tmp_User) = ResOPT;
        Res_Total_Energy = GetTotalEnergy(Num_User, Num_UAV, ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Given_P_un, Given_Q_mn_x, Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
        all_Record_VaryingK_PeakPower_Energy(current_times, tmp_User) = Res_Total_Energy;
        fprintf('PeakPower Alg. Res: %f %f(10^7); Total_Energy: %f\n', ResOPT, ResOPT/1e7, Res_Total_Energy);

        
        %% Maxmin OffloadingBits Alg
        disp('Maxmin OffloadingBits Alg');
        % Given Value
        Given_TAU_umn = ones(N, Num_User * Num_UAV) / Num_User;
        Given_L_un = ones(N, Num_User) * min(1.5 / N, CPUFreq_User*Delta /(Task_Bit_Vec(1)*c_u));
        Given_P_un = ones(N, Num_User).* Pu_max;
        Given_F_umn = ones(N, Num_User * Num_UAV) * CPUFreq_UAV / Num_User;
        Given_Q_mn_x = zeros(N,Num_UAV);
        Given_Q_mn_y = zeros(N,Num_UAV);
        Given_Qinit_mn_x = zeros(1,Num_UAV);
        Given_Qinit_mn_y = zeros(1,Num_UAV);
        init_center_loc = [MAX_X*0.25, MAX_Y*0.75; MAX_X*0.25, MAX_Y*0.25; MAX_X*0.75, MAX_Y*0.75; MAX_X*0.75, MAX_Y*0.25];
        init_r = 250;
        for m=1:Num_UAV
            Given_Qinit_mn_x(1, m) = init_center_loc(m,1) + cos(0)*init_r;
            Given_Qinit_mn_y(1, m) = init_center_loc(m,2) + sin(0)*init_r;
            for i=1:N
                tmp_theta = 2 * pi * (i/(N+1));
                Given_Q_mn_x(i, m) = init_center_loc(m,1) + cos(tmp_theta)*init_r;
                Given_Q_mn_y(i, m) = init_center_loc(m,2) + sin(tmp_theta)*init_r;
            end
        end
        % conduct alg
        compareAlg_maxminOffBits;
        if A_breakiteration ~= MAX_Iteration
            fprintf('Break! compareAlg_maxminOffBits [K:%d]-[%d]', Num_User, current_times);
            label_break = true;
            break;
        end
        % Record Res
        all_Record_VaryingK_MaxminOffBits(current_times, tmp_User) =  Record_Res_iteration(end);
        all_Record_VaryingK_MaxminOffBits_real(current_times, tmp_User) = Record_Res_real_iteration(end);
        sprintf('maxminOffBits ALG CAL_RES: %f %f', Record_Res_iteration(end), Record_Res_real_iteration(end));
        % Math Calculation
        ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
        [TargetStatic,Delay_UtilityStatic,real_Delay_UtilityStatic,prop_offloadStatic] = GetTargetValue(ck_Rate * Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, N, Num_User, Num_UAV);
        % currate Time with higher N
        fprintf('Static Math Calculation:%f', TargetStatic);
        all_Record_VaryingK_MaxminOffBits_Math(current_times, tmp_User) = TargetStatic;
        % Simulate 
        assert(N == T/Delta);
        Task_Bit_Vec;
        offload_bits = ck_Rate * Bandwidth .* Given_TAU_umn * Delta * (Matrix_Replicate_10_40');
        % diag(1:1:N) * Delta * offload_bits
        local_bits = Given_L_un * diag(Task_Bit_Vec);
        % [local_Freq]
        total_bits = offload_bits + local_bits;
        sum(total_bits) >= Task_Bit_Vec;
        tmp_res_opt = zeros(1, Num_User);
        for i = 1:N
            tmp_res_opt = tmp_res_opt + total_bits(i,:) * i;
%             tmp_res_opt = tmp_res_opt + total_bits(i,:) * i * Delta;
        end
        tmp_res_opt
        tmp_res_opt./Task_Bit_Vec;
        ResOPT = max(tmp_res_opt);
        all_Record_VaryingK_MaxminOffBits_Sim(current_times, tmp_User) = ResOPT;
        Res_Total_Energy = GetTotalEnergy(Num_User, Num_UAV, ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Given_P_un, Given_Q_mn_x, Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
        all_Record_VaryingK_MaxminOffBits_Energy(current_times, tmp_User) = Res_Total_Energy;
        fprintf('maxminOffBits Res: %f %f(10^7); Total_Energy: %f\n', ResOPT, ResOPT/1e7, Res_Total_Energy);


		%% Min Sum of Op Time Alg
        disp('Min OP Time Alg');
        % Given Value
        Given_TAU_umn = ones(N, Num_User * Num_UAV) / Num_User;
        Given_L_un = ones(N, Num_User) * min(1.5 / N, CPUFreq_User*Delta /(Task_Bit_Vec(1)*c_u));
        Given_P_un = ones(N, Num_User).* Pu_max;
        Given_F_umn = ones(N, Num_User * Num_UAV) * CPUFreq_UAV / Num_User;
        Given_Q_mn_x = zeros(N,Num_UAV);
        Given_Q_mn_y = zeros(N,Num_UAV);
        Given_Qinit_mn_x = zeros(1,Num_UAV);
        Given_Qinit_mn_y = zeros(1,Num_UAV);
        init_center_loc = [MAX_X*0.25, MAX_Y*0.75; MAX_X*0.25, MAX_Y*0.25; MAX_X*0.75, MAX_Y*0.75; MAX_X*0.75, MAX_Y*0.25];
        init_r = 250;
        for m=1:Num_UAV
            Given_Qinit_mn_x(1, m) = init_center_loc(m,1) + cos(0)*init_r;
            Given_Qinit_mn_y(1, m) = init_center_loc(m,2) + sin(0)*init_r;
            for i=1:N
                tmp_theta = 2 * pi * (i/(N+1));
                Given_Q_mn_x(i, m) = init_center_loc(m,1) + cos(tmp_theta)*init_r;
                Given_Q_mn_y(i, m) = init_center_loc(m,2) + sin(tmp_theta)*init_r;
            end
        end
        % conduct alg
        compareAlg_minOPTime;
        if A_breakiteration ~= MAX_Iteration
            fprintf('Break! compareAlg_minOPTime [K:%d]-[%d]', Num_User, current_times);
            label_break = true;
            break;
        end
        % Record Res
        % keep the optimization result
        all_Record_VaryingK_MinOPTime(current_times, tmp_User) =  Record_Res_iteration(end);
        all_Record_VaryingK_MinOPTime_real(current_times, tmp_User) = Record_Res_real_iteration(end);
        sprintf('minOPTime ALG CAL_RES: %f %f', Record_Res_iteration(end), Record_Res_real_iteration(end));
        % Math Calculation
        ck_Rate = GetAccurateRate(Given_Q_mn_x, Given_Q_mn_y, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV);
        [TargetStatic,Delay_UtilityStatic,real_Delay_UtilityStatic,prop_offloadStatic] = GetTargetValue(ck_Rate * Bandwidth, Given_TAU_umn, Given_L_un, Task_Bit_Vec, Delta, N, Num_User, Num_UAV);
        % currate Time with higher N
        fprintf('Static Math Calculation:%f', TargetStatic);
        all_Record_VaryingK_MinOPTime_Math(current_times, tmp_User) = TargetStatic;
        % Simulate 
        assert(N == T/Delta);
        Task_Bit_Vec;
        offload_bits = ck_Rate * Bandwidth .* Given_TAU_umn * Delta * (Matrix_Replicate_10_40');
        % diag(1:1:N) * Delta * offload_bits
        local_bits = Given_L_un * diag(Task_Bit_Vec);
        % [local_Freq]
        total_bits = offload_bits + local_bits;
        sum(total_bits) >= Task_Bit_Vec;
        tmp_res_opt = zeros(1, Num_User);
        for i = 1:N
            tmp_res_opt = tmp_res_opt + total_bits(i,:) * i;
%             tmp_res_opt = tmp_res_opt + total_bits(i,:) * i * Delta;
        end
        tmp_res_opt
        tmp_res_opt./Task_Bit_Vec;
        ResOPT = max(tmp_res_opt);
        all_Record_VaryingK_MinOPTime_Sim(current_times, tmp_User) = ResOPT;
        Res_Total_Energy = GetTotalEnergy(Num_User, Num_UAV, ck_Rate*Bandwidth, Given_TAU_umn, Given_L_un, Given_P_un, Given_Q_mn_x, Given_Q_mn_y, Given_Qinit_mn_x, Given_Qinit_mn_y, Task_Bit_Vec);
        all_Record_VaryingK_MinOPTime_Energy(current_times, tmp_User) = Res_Total_Energy;
        fprintf('MinSumOPTime Res: %f %f(10^7); Total_Energy: %f\n', ResOPT, ResOPT/1e7, Res_Total_Energy);
		
    end
end

save('111.mat');

filename = sprintf('VaryingKallAlg-%d-%d-%dtimes-50Mb-1500m-250m', Num_init_K, Num_end_K, Times);
tmp_filename_mat = strcat(filename, '.mat');
save(tmp_filename_mat);
% save(tmp_filename_mat,"list_num_user","all_Record_VaryingK_OPT", "all_Record_VaryingK_real_OPT","all_Record_VaryingK_Math_OPT", ...
%     "all_Record_VaryingK_Sim_OPT","all_Record_VaryingK_Sim_Local", "all_Record_VaryingK_Static", ...
%     "all_Record_VaryingK_real_Static","all_Record_VaryingK_Math_Static","all_Record_VaryingK_Sim_Static");

time_end = clock;
toc;
running_time = etime(time_end, time_begin);
fprintf('Execution Time: %f [%d-%d][%d]\n', running_time, Num_init_K, Num_end_K, Times);

y1 = mean(all_Record_VaryingK_OPT_Sim,1)*Delta;
y2 = mean(all_Record_VaryingK_Local_Sim,1)*Delta;
y3 = mean(all_Record_VaryingK_PeakPower_Sim,1)*Delta;
y4 = mean(all_Record_VaryingK_Static_Sim,1)*Delta;
y5 = mean(all_Record_VaryingK_UpperBound,1)*Delta;
y6 = mean(all_Record_VaryingK_LowerBound,1)*Delta;
y7 = mean(all_Record_VaryingK_MaxminOffBits_Sim,1)*Delta;
y8 = mean(all_Record_VaryingK_MinOPTime_Sim,1)*Delta;
%Figure

f5 = figure(5);
% plot(Num_init_K:Num_end_K, y1,'Marker','+','Color','b');
plot(list_num_user, y1,'Marker','o','LineWidth',2,'Color','b');
hold on;
plot(list_num_user, y2,'Marker','+','LineWidth',2,'Color','k');
plot(list_num_user, y3,'Marker','*','LineWidth',2,'Color','r');
plot(list_num_user, y4,'Marker','d','LineWidth',2,'Color','g');
plot(list_num_user, y5, '--', 'Marker','^','LineWidth',2,'Color','g');
plot(list_num_user, y6, '--', 'Marker','v','LineWidth',2,'Color','k');
plot(list_num_user, y7, '--', 'Marker','d','LineWidth',2,'Color','y');
plot(list_num_user, y8, '-.', 'Marker','*','LineWidth',2,'Color','r');
xlabel('K');
% ylabel('Min-max AoT');
ylabel('Max AoT of GUs');
legend('AoT OPT Alg.','Local Comp.', 'Peak Power', 'Static UAV', 'Upper Bound', 'Lower Bound', 'Maxmin OffBits Alg.', 'Min OP Time Alg.');
tmp_filename_jpg = strcat(filename, '.jpg');
saveas(f5, tmp_filename_jpg);
% tmp_filename_eps = strcat(filename, '.eps');
% saveas(f5, tmp_filename_eps);


y1 = mean(all_Record_VaryingK_OPT_Energy,1);
y2 = mean(all_Record_VaryingK_Local_Energy,1);
y3 = mean(all_Record_VaryingK_Static_Energy,1);
y4 = mean(all_Record_VaryingK_PeakPower_Energy,1);
y5 = mean(all_Record_VaryingK_MaxminOffBits_Energy,1);
y6 = mean(all_Record_VaryingK_MaxminOffBits_Energy,1);
f6 = figure(6);
% clf;
hold on;
idx_themes = 0;
pos_idxC = idx_themes * 6;
p = plot(list_num_user, y1,'-','Marker','o','LineWidth',2,'Color','b');
% set(p, 'markerfacecolor', get(p, 'color'));
p = plot(list_num_user, y2,'-','Marker','s','LineWidth',2,'Color','r');
% set(p, 'markerfacecolor', get(p, 'color'));
p = plot(list_num_user, y3,'-','Marker','>','LineWidth',2,'Color','k');
% set(p, 'markerfacecolor', get(p, 'color'));
p = plot(list_num_user, y4,':','Marker','+','LineWidth',2,'Color','g');
% set(p, 'markerfacecolor', get(p, 'color'));
p = plot(list_num_user, y5,'--','Marker','d','LineWidth',2,'Color','r');
% set(p, 'markerfacecolor', get(p, 'color'));
p = plot(list_num_user, y6,'--','Marker','hexagram','LineWidth',2,'Color','b');
% set(p, 'markerfacecolor', get(p, 'color'));
% xlabel('K');
ylim([0, 120000]);
xlabel('Number of GUs');
ylabel('Total Energy Consumption (J)');
l = legend('AoT OPT Alg.','Local Comp.', 'Static UAV', 'Peak Power', 'MaxOffBits', 'MinOPtime');
%l.Position = [0.63,0.74,0.1,0.1];
l.NumColumns = 3;
l.Location = 'North';
grid on;
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
box on;

diary off;
fprintf('Now:');
fprintf(' %d',clock);
fprintf('\n');