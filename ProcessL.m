function Given_L_un = ProcessL(ck_Rate, Given_TAU_umn, Given_L_un, Task_Bit_Vec, delta, CPUFreq_User, c_u, N, Num_User, Num_UAV, Bandwidth)
%ProcessL Get the accurate Target Value
%   此处显示详细说明
Matrix_Replicate_10_40 = zeros(Num_User, Num_User * Num_UAV);
for u=1:Num_User
    for m=1:Num_UAV
        Matrix_Replicate_10_40(u, u+(m-1)*Num_User) = 1;
    end
end

ck_offloadbits = sum((ck_Rate .* Given_TAU_umn) * delta * Matrix_Replicate_10_40', 1);
ck_localbits = sum(Given_L_un * diag(Task_Bit_Vec/Bandwidth),1);
ck_bits = ck_offloadbits + ck_localbits; 
for user_i=1:Num_User
    if ck_bits(user_i) < (Task_Bit_Vec(user_i)/Bandwidth)
        left_prop = ( Task_Bit_Vec(user_i)/Bandwidth-ck_bits(user_i) )/(Task_Bit_Vec(user_i)/Bandwidth);
        assert (left_prop > 0)
        max_capacity = (CPUFreq_User * delta)/(Task_Bit_Vec(user_i)*c_u);
        for time_j=1:N
            if left_prop <= 0
                break
            else
                if Given_L_un(time_j, user_i) < max_capacity
                    old_value = Given_L_un(time_j, user_i);
                    Given_L_un(time_j, user_i) = min(max_capacity, old_value + left_prop);
                    increment = Given_L_un(time_j, user_i) - old_value;
                    left_prop = left_prop - increment;
                end
            end
        end
    end
end
end

