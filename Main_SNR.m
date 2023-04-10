clear;
%% parameters of methods- ista,omp,sbl
n_batch = 100; % number of batches
lambda = 5; % regulation ratio of optimazation problem
max_iter_ista = 1000; % maximum iterations of ista
sparsity_level_omp = 100; % maximum iterations of omp
max_iter_sbl = 100; % maximum iterations of sbl
%% variable area
% how many group of signals we observe
% num_signals = [5 10 20 30 40];
for num_signals = 12
    SNR = [-10 -5 0 5 10 20]; % SNR value in dB
    recorder_mean = zeros(4,length(SNR));
    recorder_varn = zeros(4,length(SNR));
    recorder_time = zeros(4,length(SNR));
    % assume signal power is 1, noise power is 1/SNR
    sigma = sqrt(1./(10.^(SNR./10)));
    G=40;
    %% mean error test for three methods
    for j=1:length(SNR)
        error_list_ista = zeros(n_batch, 1); % list of error for each batch
        error_list_omp = zeros(n_batch, 1);
        error_list_sbl = zeros(n_batch, 1);
        error_list_sblu = zeros(n_batch, 1);

        % the mean of ista
        tStart1 = cputime;
        for i = 1:n_batch
            % generate channel
            [y_bar,H,Q_bar,noise,A_R_bar,A_T_bar,H_vec] = channel_generate(sigma(1,j), num_signals);
            % carry out ista
            H_v = ista_complex(Q_bar,y_bar, lambda, max_iter_ista);
            % calculate H_est
            H_a_est_full = reshape(H_v, [G, G]);
            H_est = A_R_bar * H_a_est_full*A_T_bar';
            % insert error into error_list
            error_ista = norm(H_est-H)/norm(H);
            error_list_ista(i) = error_ista^2;
        end
        tEnd1 = cputime - tStart1;
        % mean
        disp("mean of ISTA: "+string(mean(error_list_ista)))
        % standard deviation
        disp("standard deviation of ISTA: "+string(std(error_list_ista)))
        recorder_mean(1,j) = mean(error_list_ista);
        recorder_varn(1,j) = std(error_list_ista);
        recorder_time(1,j) = tEnd1/n_batch;
        % the mean of omp
        tStart2 = cputime;
        for i = 1:n_batch
            % generate channel
            [y_bar,H,Q_bar,noise,A_R_bar,A_T_bar,H_vec] = channel_generate(sigma(1,j), num_signals);
            % carry out omp
            H_v = omp_complex(Q_bar,y_bar,sparsity_level_omp);
            % calculate H_est
            H_a_est_full = reshape(H_v, [G, G]);
            H_est = A_R_bar * H_a_est_full*A_T_bar';
            % insert error into error_list
            error_omp = norm(H_est-H)/norm(H);
            error_list_omp(i) = error_omp^2;
        end
        tEnd2 = cputime - tStart2;
        % mean
        disp("mean of OMP: "+string(mean(error_list_omp)))
        % standard deviation
        disp("standard deviation of OMP: "+string(std(error_list_omp)))
        recorder_mean(2,j) = mean(error_list_omp);
        recorder_varn(2,j) = std(error_list_omp);
        recorder_time(2,j) = tEnd2/n_batch;
        % the mean of sbl
        tStart3 = cputime;
        for i = 1:n_batch
            % generate channel
            [y_bar,H,Q_bar,noise,A_R_bar,A_T_bar,H_vec] = channel_generate(sigma(1,j), num_signals);
            % carry out sbl
            H_v = SBL(y_bar,Q_bar,sigma(1,j),1e-4,max_iter_sbl);
            % calculate H_est
            H_a_est_full = reshape(H_v, [G, G]);
            H_est = A_R_bar * H_a_est_full*A_T_bar';
            % insert error into error_list
            error_sbl = norm(H_est-H)/norm(H);
            error_list_sbl(i) = error_sbl^2;
        end
        tEnd3 = cputime - tStart3;
        % mean
        disp("mean of SBL: "+string(mean(error_list_sbl)))
        % standard deviation
        disp("standard deviation of SBL: "+string(std(error_list_sbl)))
        recorder_mean(3,j) = mean(error_list_sbl);
        recorder_varn(3,j) = std(error_list_sbl);
        recorder_time(3,j) = tEnd3/n_batch;
        % the mean of sbl (uknown sigma)
        tStart4 = cputime;
        for i = 1:n_batch
            % generate channel
            [y_bar,H,Q_bar,noise,A_R_bar,A_T_bar,H_vec] = channel_generate(sigma(1,j), num_signals);
            % carry out sblu
            H_v = SBLU(y_bar,Q_bar,noise,1e-4,max_iter_sbl);
            % calculate H_est
            H_a_est_full = reshape(H_v, [G, G]);
            H_est = A_R_bar * H_a_est_full*A_T_bar';
            % insert error into error_list
            error_sblu = norm(H_est-H)/norm(H);
            error_list_sblu(i) = error_sblu^2;
        end
        tEnd4 = cputime - tStart4;
        % mean
        disp("mean of SBLU: "+string(mean(error_list_sblu)))
        % standard deviation
        disp("standard deviation of SBLU: "+string(std(error_list_sblu)))
        recorder_mean(4,j) = mean(error_list_sblu);
        recorder_varn(4,j) = std(error_list_sblu);
        recorder_time(4,j) = tEnd4/n_batch;
    end
    % plot the mean CPU time
    figure(1);
    bar(SNR,10*log10(recorder_time'));
    title(['Mean CPU Time # pilot = ',num2str(num_signals)]);
    xlabel("SNR/dB");
    ylabel("CPU Time/dBs");
    legend('ISTA','OMP','SBL','SBLU');
    saveas(gcf, ['CPUTime',num2str(num_signals),'pilots.png']);
    % plot the normalized mean square error
    figure(2);
    plot(SNR,recorder_mean(1,:),'-o',SNR,recorder_mean(2,:),'-o',SNR,recorder_mean(3,:),'-o',SNR,recorder_mean(4,:),'-o');
    title(['Normalized Mean Square Error # pilot = ',num2str(num_signals)]);
    xlabel("SNR/dB");
    ylabel("MSE");
    legend('ISTA','OMP','SBL','SBLU');
    saveas(gcf, ['MSE',num2str(num_signals),'pilots.png']);
    % plot the variance of normalized mean square error
    figure(3);
    plot(SNR,recorder_varn(1,:),'-o',SNR,recorder_varn(2,:),'-o',SNR,recorder_varn(3,:),'-o',SNR,recorder_varn(4,:),'-o');
    title(['Variance of MSE # pilot = ',num2str(num_signals)]);
    xlabel("SNR/dB");
    ylabel("Variance");
    legend('ISTA','OMP','SBL','SBLU');
    saveas(gcf, ['Variance',num2str(num_signals),'pilots.png']);
end