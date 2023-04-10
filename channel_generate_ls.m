function [y,H,F] = channel_generate_ls(sigma,num_obs)
% this version is only used for Least Square. 
% y is the recived/observation matrix, H is the channel matrix, F is the measurement matrix. 
% MIMO channel estimation
% Here is a code tells you how to generate the channel matrix. You need to
% estimate the channel with your own proposed techniques.
% The prerequisite for understanding this code is understanding the channel
% part of the paper. If anything is unclear with channel model, feel free
% to contact me.

% Some pre-defined values
N_R = 32; % the number of antennas in BS
N_T = 16; % the number of antennas in User
N_RF = 4; % the number of RF links

% Design the dictionary
N = 40; % the number of grids

idx_t=randsample([1:N],N_RF);
idx_r=randsample([1:N],N_RF);


% generate AoDs Angles of Departure

AoD_ms = (idx_t-1)*pi/N;
% generate AoAs Angles of Arrival
AoA_bs = (idx_r-1)*pi/N;
% generate path gains CN(0,1)
g = (1 * randn(N_RF, 1) + 1i * randn(N_RF, 1)) / sqrt(2);
% generate channel
H = zeros(N_R, N_T);

for i = 1:N_RF
    H = H + g(i) * generate_steering(N_R, AoA_bs(i)) * generate_steering(N_T, AoD_ms(i))';
end

H = H * sqrt(N_R * N_T / N_RF);
%%
% some assumptions written here
% transmit power is 1
P=1;

% how many signals are transmitted
N_beam = num_obs;

% core formula
% y = sqrt(P) ((A_T' * F)' \otimes A_R) * vec(H) + noise
% noise power follows complex Gaussian with mu=0 and var=sigma^2

N_vec = N_RF * N_RF;
% generate noise, N_R*N_beam entries
noise = (1 * randn(N_R * N_beam, 1) + 1i * randn(N_R * N_beam, 1)) / sqrt(2) * sigma;
% generate transmit signal
% the signal has N_T entries and 5 columns, s is a matrix, gaussian
F = (1 * randn(N_T, N_beam) + 1i * randn(N_T, N_beam)) / sqrt(2);


A_T = zeros(N_T, N_RF);
for i = 1:N_RF
    A_T(:,i) = [generate_steering(N_T, AoD_ms(i))];
end

A_R = zeros(N_R, N_RF);
for i = 1:N_RF
    A_R(:,i) = [generate_steering(N_R, AoA_bs(i))];
end


G = N;
A_T_bar = zeros(N_T, G);
for i = 1:G
    A_T_bar(:,i) = [generate_steering(N_T, (i-1)*pi/G)];
end

A_R_bar = zeros(N_R, G);
for i = 1:G
    A_R_bar(:,i) = [generate_steering(N_R, (i-1)*pi/G)];
end

H_a = sqrt(N_R * N_T / N_RF)*diag(g);

% % calculate (A_T' * F)'
A_T_F = (A_T_bar' * F).';

% calculate Q = (A_T' * F)' \otimes A_R
Q_bar = kron(A_T_F, A_R_bar);

% diagonal matrix of size N_R*N_R
W = diag(ones(N_R, 1));
Q = kron(F.',W);
H_vec = H(:);

% generate y
Y = H*F;
y_bar = Y(:)+noise;
y = reshape(y_bar, [32, num_obs]);

%% function space
function array_res = generate_steering(dim, angle)
    % generate the array response of the antenna array
    % dim: dimension of the vector
    % angle: AoA or AoD
    array_res = zeros(dim, 1);

    for i_ = 1:dim
        array_res(i_, :) = sqrt(1 / dim) * exp(1j * 2 * pi * 0.5 * (i_ - 1) * cos(angle));
    end

end

end