%% Calculate log_sensitivity for rings via kde (calc_log_sens_ring_kde.m)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function takes an input the controller data for spin rings and dephasing operators and
% calls on the AnalyzeRobustnessDephasing.m routine to compute the kde-based log-sensitivity.  

% On input:
% N - ring size 
% out - target spin
% opt - optimization option (dephasing, fidelity, or overlap)

% On output - saved to ../results/kde-rings/log_kde_opt_N-ring_1-out.mat
% density - density of the fidelity error distribution computed by the AnalyzeErrorDensity.m routine   
% Error - array of error values computed by the AnalyzeRobustnessDephasing.m

close all; clearvars -except noise_choice noise

if ~exist('../results/kde-rings/','dir')
    mkdir('../results/kde-rings/');
end

num=1000;
in = 1;
option = {'dephasing';'fidelity';'overlap'};

for x = 1:3
    opt = option{x};
	for N = 5:6
	    dop = sprintf('../results/dephasing_op/dephasingop_ld_ring-xx-%d',N);
	    load(dop);
	    deph = arrayfun(@(n) DephasingOp{n},1:num,'UniformOutput',false);

	    for out = 2:floor(N/2)+1     

			if exist(sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d.mat',opt,N,out)) == 2
				disp(sprintf('../results/kde-rings/log_kde_%s_%d_ring_1-%d exists, skipping',opt,N,out))
			else
		        disp(sprintf('Optimzer=%s N=%d out=%d %s progress:',opt,N,out))
		        loadtag = sprintf('../results/controllers-rings/%s_%d-ring_1-%d.mat',opt,N,out);
		        load(loadtag);
     
		        % Inner loop for computation of data for each controller
		        for k = 1:100
		            status = sprintf(' => Optimizer=%s N=%d target=%d controller # = %d',opt,N,out,k);
		            disp(status) 
		            controller = sys{k};

		            % Inner loop for computation per perturbation  
		            for run = 1:num   
		                error{run} = AnalyzeRobustnessDephasing(controller,deph{run},in,out);           
		            end
		            Error{k} = cell2mat(error');
		            density{k} = AnalyzeErrorDensity(Error{k},controller);
		            close all;
		        end
                
		        % save results 
		        savetag = sprintf('../results/kde-rings/log_kde_%s_%d-ring_1-%d',opt,N,out);
		        save(savetag,'density','Error');
		    end
		end
	end
end

