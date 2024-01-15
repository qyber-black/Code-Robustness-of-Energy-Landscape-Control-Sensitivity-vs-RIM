%% Calculate log_sensitivity for chains analytically (calc_log_sens_chain.m)

function calc_log_sens_chain(noise_level)
%
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 S M Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 
%
% This function takes as input the controller data for spin chains and 
% dephasing operators and calculates the sensitivity and log-sensitivity 
% for dephasing analytically. 
% 
% Input: noise level from the following list of options
% {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
% default if not specified '0.1'
% 
% On input:
% N - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

% On output     - saved to ../results/log_sens-chains/log_sens_0.1_opt_N-chain_out.mat
% mean_log_sens - mean log-sensitivity across all dephasing operators 
% log_sens      - total log-sensitivity array 
% sens          - array of sensitivity values 
% fid_rim       - fidelity based on RIM data
% fid_calc      - fidelity calculated analytically 
% err           - fidelity error from fid_calc

% select noise level (default 0.1)
noise={'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};

if ~exist('noise_level','var') || all(cellfun(@(x) strcmp(x,noise_level), noise))
    noise_level = '0.1'
end

disp(sprintf('noise level selected %s',noise_level))

if strcmp(noise_level,"0.0")
    q0 = 1; % include l-bfgs-controllers
else
    q0 = 2; % no l-bfgs controllers
end

if ~exist('../results/log_sens-chains','dir')
    mkdir('../results/log_sens-chains');
end

num=1000;

in = 1;

type = {'lbfgs','nmplus','ppo','snob'};

for N = 5:6

    b = bloch_basis(N,eye(N)); % produce basis matrices for density matrix expansion
    off = (N^2-N)/2;

    % load dephasing operators
    dop = sprintf('../results/dephasing_op/dephasingop_ld_ring-xx-%d.mat',N);
    load(dop);
    deph = arrayfun(@(n) DephasingOp{n},1:num,'UniformOutput',false);

    % compute S_mu based on dephasing operators 
    for ell = 1:1000
        G = deph{ell}+deph{ell}';
        G = G/sum(G(1:end));
        S{ell} = zeros(N^2,N^2);
        count = 1;
        for m = 1:N-1
            for n = m+1:N
                S{ell}(count,count) = G(m,n);
                S{ell}(count+off,count+off) = G(m,n);
                count = count+1;
            end
        end
    end

    if N == 5
        target = [3 5];
    else
        target = [4 6];
    end
    for x = 1:2
        out = target(x);
        for q = q0:4
            opt = type{q};
            % Check for existence of data and skip if already present  
            if exist(sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out)) == 2
               disp(sprintf('Noise=%s Optimizer=%s N=%d out=%d complete',noise_level,opt,N,out));
            else 
               disp(sprintf('Noise=%s Optimizer=%s N=%d out=%d %s progress:',noise_level,opt,N,out))
               loadtag = sprintf('../results/controllers-chains/noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
               load(loadtag);

               X = zeros(N^2,N^2); 

               % Inner loop for computation of data for each controller
               for k = 1:length(sys)
                   status = sprintf('Noise=%s Optimizer=%s N=%d target=%d controller # = %d',noise_level,opt,N,out,k);
                   disp(status) 
                   D  = diag(sys{k}.x); % produce control matrix
                   Hd = sys{k}.H;       % produce controlled Hamiltonian 
                   t  = sys{k}.T;       % extract read-out time 
                   r0 = zeros(N); r0(in,in)  =1;   % initial state 
                   r1 = zeros(N); r1(out,out)=1;   % target state 
                   [A,lambda,V,r_in,r_out] = bloch(N,Hd,r0,r1,b); 
                   fid_rim(k,1)  = sys{k}.fidelity;
                   fid_calc(k,1) = r_out'*expm(A*t)*r_in;
                   Diff(k,1) = fid_rim(k,1) - fid_calc(k,1); 
                   err = 1-fid_rim;

                   % Inner loop for computation per perturbation  
                   for run = 1:num        
                       X = S{run}*t*expm(A*t);
                       xi = 1;
                       sens(k,run) = -r_out'*X*r_in;  % compute sensitivity 
                       log_sens(k,run) = abs(xi*(-1)*r_out'*X*r_in)/(err(k,1)); % compute log-sensitivity            
                   end

                   log_sens(k,num+1) = mean(log_sens(k,1:num)); % compute norm of log-sensitivity 
                   sens(k,num+1) = mean(sens(k,1:num));
                   mean_log_sens = log_sens(:,1001);
               end

               % save results 
               savetag = sprintf('../results/log_sens-chains/log_sens_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out)
               save(savetag,'mean_log_sens','log_sens','sens','fid_rim','fid_calc','Diff','err');
            end
        end
    end
end
