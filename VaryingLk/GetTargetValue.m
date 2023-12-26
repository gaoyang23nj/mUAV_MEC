function [Target,Delay_Utility,real_Delay_Utility,prop_offload] = GetTargetValue(ck_Rate, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, N, Num_User, Num_UAV)
% GetTargetValue  get the target value
Matrix_Replicate_10_40 = zeros(Num_User, Num_User * Num_UAV);
for u=1:Num_User
    for m=1:Num_UAV
        Matrix_Replicate_10_40(u, u+(m-1)*Num_User) = 1;
    end
end

offload_bits = ck_Rate/1e6 .* Given_TAU_umn * delta * (Matrix_Replicate_10_40');
offload_ratio = offload_bits ./ repmat(Task_Bit_Vec/1e6,[N,1]);
Delay_Utility = sum((diag(1:1:N)) * (offload_ratio + Given_L_un), 1);
Target = max(Delay_Utility);

prop_offload = ones(1,N) * offload_ratio;
prop_local = ones(1,N) * Given_L_un;

real_offload_bits = ck_Rate .* Given_TAU_umn * delta * (Matrix_Replicate_10_40');
real_local_bits = Given_L_un * diag(Task_Bit_Vec);
real_Delay_Utility = sum((diag(1:1:N)) * (real_offload_bits + real_local_bits), 1);

fprintf('Target:%d\n', Target);
% Delay_Utility;
% real_Delay_Utility;
end

