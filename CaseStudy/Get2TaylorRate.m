function ck_Rate_2taylor_true = Get2TaylorRate(Given_Q_mn_x, Given_Q_mn_y, Given_Q_mn_x_r, Given_Q_mn_y_r, Loc_User_x, Loc_User_y, Given_P_un, H, Sigma2, rho, N, Num_User, Num_UAV)
%GetAccurateRate get the accurate Communication Rate (without Bandwidth), for all u,m,n
%  log

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

ck_uav_x = Given_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
ck_uav_y = Given_Q_mn_y * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
ck_denomin = power(ck_uav_x, 2) + power(ck_uav_y, 2) + H*H;
% Rate hat
ck_Rate_hat_lg = ((Given_P_un * rho * Matrix_Replicate_10_40) ./ ck_denomin) * Matrix_Replicate_4_40' + Sigma2;
ck_Rate_hat = log2(ck_Rate_hat_lg);
% Rate tilde
ck_Rate_tilde_lg = ((Given_P_un .* rho * Matrix_Replicate_10_40) ./ ck_denomin) * Matrix_delete_user + Sigma2;
ck_Rate_tilde = log2(ck_Rate_tilde_lg);
ck_Rate = ck_Rate_hat * Matrix_Replicate_4_40 - ck_Rate_tilde;

dist_uav_user_x_r = Given_Q_mn_x_r * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
dist_uav_user_y_r = Given_Q_mn_y_r * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
denomin_r = power(dist_uav_user_x_r, 2) + power(dist_uav_user_y_r, 2) + H*H;
% Rate hat Taylor
Rate_hat_Taylor_lg = (((Given_P_un * rho * Matrix_Replicate_10_40 / Sigma2) ./ denomin_r  ) * Matrix_Replicate_4_40') + 1;
Rate_hat_Taylor_term1 = log2(Rate_hat_Taylor_lg) + log2(Sigma2);
Rate_hat_Taylor_mul = ((Given_P_un * rho * Matrix_Replicate_10_40 * log2(exp(1)) / Sigma2) ./ power(denomin_r, 2)) ./ (Rate_hat_Taylor_lg * Matrix_Replicate_4_40);
%dist_uav_user_x = Var_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_x * Matrix_Replicate_10_40;
%dist_uav_user_y = Var_Q_mn_x * Matrix_Replicate_4_40 - ones(N,1) * Loc_User_y * Matrix_Replicate_10_40;
dist_uav_uav_x_r = Given_Q_mn_x * Matrix_Replicate_4_40 - Given_Q_mn_x_r * Matrix_Replicate_4_40;
dist_uav_uav_y_r = Given_Q_mn_y * Matrix_Replicate_4_40 - Given_Q_mn_y_r * Matrix_Replicate_4_40;
Rate_hat_Taylor_mul_dist = ((dist_uav_user_x_r .* dist_uav_uav_x_r) + (dist_uav_user_y_r .* dist_uav_uav_y_r))*2;
Rate_hat_Taylor_term2 = (Rate_hat_Taylor_mul .* Rate_hat_Taylor_mul_dist) * Matrix_Replicate_4_40';
Rate_hat_Taylor = Rate_hat_Taylor_term1 - Rate_hat_Taylor_term2;
% Rate tilde Taylor
Rate_tilde_Taylor_lg =  (((Given_P_un * rho * Matrix_Replicate_10_40 / Sigma2) ./ denomin_r) * Matrix_delete_user) + 1;
Rate_tilde_Taylor_term1 = log2(Rate_tilde_Taylor_lg) + log2(Sigma2);
Rate_tilde_Taylor_mul = ((Given_P_un * rho * Matrix_Replicate_10_40 * log2(exp(1)) / Sigma2) ./ power(denomin_r, 2)) ./ Rate_tilde_Taylor_lg;
Rate_tilde_Taylor_term2 = (Rate_tilde_Taylor_mul .* Rate_hat_Taylor_mul_dist) * Matrix_delete_user;
Rate_tilde_Taylor = Rate_tilde_Taylor_term1 - Rate_tilde_Taylor_term2;
            
ck_Rate_2taylor_true = Rate_hat_Taylor * Matrix_Replicate_4_40 - Rate_tilde_Taylor;
% ck_Rate = ck_Rate * Bandwidth;

A_rate_error = power(ck_Rate - ck_Rate_2taylor_true, 2);
end

