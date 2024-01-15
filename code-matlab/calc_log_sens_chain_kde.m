%% Calculate log_sensitivity for chains via kde (calc_log_sens_chain_kde.m)

function [] = calc_log_sens_chain_kde(noise_level)

% SPDX-License-Identifier: CC-BY-SA-4.0 
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 S M Shermer <lw1660@gmail.com>

% This function takes an input the controller data for spin chains and dephasing operators and
% calls on the AnalyzeRobustnessDephasing.m routine to compute the kde-based log-sensitivity.  

% Input: noise level from the following list of options
% {'0.0','0.01','0.02','0.03','0.04','0.05','0.1'};
% default if not specified '0.1'

% On input:
% N   - chain size 
% out - target spin
% opt - optimization option (nmplus, ppo, or snob)

% On output - saved to ../results/log_sens-chains/log_kde_0.1_opt_N-chain_out.mat
% density   - density of fidelity error distribution computed by AnalyzeErrorDensity.m   
% Error     - array of error values computed by AnalyzeRobustnessDephasing.m

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

if ~exist('../results/kde-chains/','dir')
    mkdir('../results/kde-chains/');
end

num=1000;
in = 1;
type = {'lbfgs','nmplus','ppo','snob'};

for N = 5:6
	dop = sprintf('../results/dephasing_op/dephasingop_ld_ring-xx-%d',N);    
	load(dop);
	deph = arrayfun(@(n) DephasingOp{n},1:num,'UniformOutput',false);
	
	if N == 5		
		target = [3 5];	
	else
        target = [4 6];    
	end

    for x = 1:2
		out = target(x);

	    for q = q0:length(type);
	        opt = type{q};

	        if exist(sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out)) == 2
	            disp(sprintf('Noise=%s Optimizer=%s N=%d out=%d %s complete',noise_level,opt,N,out))
	        else 
		        disp(sprintf('Noise=%s Optimizer=%s N=%d out=%d %s progress:',noise_level,opt,N,out))
	    	    %loadtag = sprintf('../data-raw/controllers-irtaza/noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
	    	    loadtag = sprintf('../results/controllers-chains/noise_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out);
                load(loadtag);

	    	    % Inner loop for computation of data for each controller
	    	    for k = 1:length(sys)
	    	        status = sprintf('Noise=%s Optimizer=%s N=%d target=%d controller # = %d',noise_level,opt,N,out,k);
	    	        disp(status) 
	    	        controller = sys{k};

	    	        % Inner loop for computation per perturbation  
	    	        for run = 1:num   
	    	            error{run} = AnalyzeRobustnessDephasing(controller,deph{run},in,out);           
	    	        end
    		        Error{k} = cell2mat(error');
    		        density{k} = AnalyzeErrorDensity(Error{k},controller);
    		    end

    		    % save results 
    		    savetag = sprintf('../results/kde-chains/log_kde_%s_%s_%d-chain_1-%d.mat',noise_level,opt,N,out)
    	    	save(savetag,'density','Error');
    	    
    		end
		end
	end
end
            


