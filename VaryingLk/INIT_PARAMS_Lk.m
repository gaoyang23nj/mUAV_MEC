%% Params
T = 100;
N = 200;
% delta = 0.5;
Delta = T/N;
% altidue 100m
H = 100;
% max velocity, 30m/s
v_max = 30;
% min distance between two UAVs, 50m
d_min = 100;
% Num_UAV = 4;
% Num_User = 10;
MAX_X = 1500;
MAX_Y = 1500;
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
E_uav_OE_max = 2000;
E_uav_prop_max = 50000;

% % !!! the random seed
% rng_seed = 1000;
% rng(rng_seed);
% t = clock;
% rng(t(6));
% rand(1)
   
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


