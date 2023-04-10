clear;
%% parameters of methods- ista,omp,sbl
n_batch = 1000; % number of batches
lambda = 5; % regulation ratio of optimazation problem
max_iter_ista = 1000; % maximum iterations of ista
sparsity_level_omp = 100; % maximum iterations of omp
max_iter_sbl = 100; % maximum iterations of sbl
%% variable area
% SNR = [-5 0 10 20]
for SNR = [10,20] 
    % SNR = 10; % SNR value in dB
    num_signals = [16 20 28]; % how many group of signals we observe
    recorder_mean = zeros(1,length(num_signals));
    recorder_varn = zeros(1,length(num_signals));
    recorder_time = zeros(1,length(num_signals));
    % assume signal power is 1, noise power is 1/SNR
    sigma = sqrt(1/(10^(SNR/10)));
    G=40;
    %% mean error test for three methods
    for j=1:length(num_signals)
        error_list_ls = zeros(n_batch, 1); % list of error for each batch

        % the mean of ls
        tStart1 = cputime;
        for i = 1:n_batch
            % generate channel
            [y,H,Q] = channel_generate_ls(sigma, num_signals(1,j));
            % carry out LS (Least Square)
            H_est= y*Q'*(Q*Q')^(-1); 
            % insert error into error_list
            error_ls = norm(H_est-H)/norm(H);
            error_list_ls(i) = error_ls^2;
        end
        tEnd1 = cputime - tStart1;
        % mean
        disp("mean of LS: "+string(mean(error_list_ls)))
        % standard deviation
        disp("standard deviation of LS: "+string(std(error_list_ls)))
        recorder_mean(1,j) = mean(error_list_ls);
        recorder_varn(1,j) = std(error_list_ls);
        recorder_time(1,j) = tEnd1/n_batch;
    end
    % plot the mean CPU time
    figure(1);
    bar(num_signals,recorder_time');
    title(['Mean CPU Time (SNR = ',num2str(SNR),' dB)']);
    xlabel("Pilot Number");
    ylabel("CPU Time/s");
    legend('LS');
    saveas(gcf, ['CPUTime',num2str(SNR),'dB.png']);
    % plot the normalized mean square error
    figure(2);
    plot(num_signals,recorder_mean(1,:),'-o');
    title(['Normalized Mean Square Error (SNR = ',num2str(SNR),' dB)']);
    xlabel("Pilot Number");
    ylabel("MSE");
    legend('LS');
    saveas(gcf, ['MSE',num2str(SNR),'dB.png']);
    % plot the variance of normalized mean square error
    figure(3);
    plot(num_signals,recorder_varn(1,:),'-o');
    title(['Variance of MSE (SNR = ',num2str(SNR),' dB)']);
    xlabel("Pilot Number");
    ylabel("Variance");
    legend('LS');
    saveas(gcf, ['Variance',num2str(SNR),'dB.png']);
end