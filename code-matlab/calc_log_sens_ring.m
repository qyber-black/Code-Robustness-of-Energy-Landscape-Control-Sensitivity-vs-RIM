%% Calcuate log-sensitivy from LTI formulation for rings
function [] = calc_log_sens_ring()

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This routine takes as input a set of 100 controllers optimized for either
% fidelity, overlap, or transfer under dephasing along with a set of 1000
% dephasing operators. It outputs a data file with the sensitivity and
% log-sensitiviy calculated from the analytical expression as well as plots
% of the log-sensitivity versus fidelity error. 

% On Input or Used for Execution Computations
% opt - optimization option (dephasing, fidelity, or overlap)
% N - size of spin-ring
% S - structure of 1000 perturbation (dephasing) matrices 
% out - target spin for transfer (1 -> out) 
% err - nominal error extracted from controllers_opt_ld_ring-xx-N_0-out.mat
%       file 
% time / t - time of excitation transfer 
% D - diagonal matrix of control bias fields 
% Hd - controlled Hamiltonian (H+D)
% v - matrix of eigenvalues of Hd (Hd = v*diag(lambda)*v')
% lambda - eigenvalues of Hd 
% omega - skew-symmetric matrix of transition frequencies 
% Omega - vecorized representation of Omega
% A - state matrix for LTI formulation (A = -i*Omega)
% r0 - initial density matrix in Hamiltonian basis
% r1 - target density matrix in Hamiltonian basis 
% r_in - vectorized representation of r0
% r_out - vectorized representation of r1 
% err_calc - nominal error calculation for LTI formulation  
% diff - difference in err and err_calc
%
% Output is saved to
% ../results/log_sens-rings/log_sens_opt_N_out:
% log_sens - array of log-sensitivity of the error
%            each row is the results for a given controller 
%            each column shows the results by dephasing operator 
%            the final column is the mean of the error across all 1000
%            dephasing operators
% sens - array of differential sensitivity ordered same as log-sensitivity 
% err  - described above 
% err_calc - described above 
% diff - described above                  

close all; clearvars -except noise_choice noise

in  = 1;
num = 1000;
tol = 10^-12;
option = {'dephasing';'fidelity';'overlap'};

if ~exist('../results/log_sens-rings','dir')
    mkdir('../results/log_sens-rings');
end

% Loop for each optimization option 
for x = 1:3
    opt = option{x};

	% Loop for each ring size
	for N = 5:6
	    off = (N^2-N)/2;
	    I = eye(N);
	    b = bloch_basis(N,I);
	    e = arrayfun(@(n) I(:,n),1:N,'UniformOutput',0);
    
	    % Build perturbation structure for each dephasing matrix 
	    tag = sprintf('../results/dephasing_op/dephasingop_ld_ring-xx-%d',N);
	    load(tag);
	    deph = arrayfun(@(n) DephasingOp{n},1:num,'UniformOutput',false);
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
   
	    % Loop for each target spin
    	for out = 2:floor(N/2)+1
    
	    	% Check for existence of data for transfer and skip if already present  
	    	if exist(sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d.mat',opt,N,out)) == 2
	        	disp(sprintf('N=%d out=%d %s complete',N,out,opt))
	    	else 
				disp(sprintf('N=%d out=%d %s progress:',N,out,opt))
				X = zeros(N^2,N^2);   % Initialize X for computation of log-sens 
				tag1 = sprintf('../results/controllers-rings/%s_%d-ring_1-%d',opt,N,out);
				load(tag1); % Load controller data files 
	       		err = arrayfun(@(x) 1-sys{x}.fidelity,1:100)'; %extract nominal error (no perturbations)
	       		time = arrayfun(@(x) sys{x}.T,1:100)'; 
  
	       		% Inner loop for computation of data for each controller
	       		M = max(size(err));
            
	            for k = 1:M
					status = sprintf('N=%d target=%d option=%s controller # = %d',N,out,opt,k);
				    disp(status)    

					D = diag(sys{k}.x);	 % produce control matrix
					Hd = sys{k}.H;   	 % produce controlled Hamiltonian 		           
					r0 = zeros(N); r0(in,in)  =1; 
					r1 = zeros(N); r1(out,out)=1; 	           
					[A,lambda,V,r_in,r_out] = bloch(N,Hd,r0,r1,b);
					t = time(k,1);				% extract transfer time	     	      
					err_calc(k,1) = 1-r_out'*expm(t*A)*r_in;
           
					% Inner loop for computation per perturbation  
					for run = 1:num        
		              X = -S{run}*t*expm(-t*A);
		              xi = 1;
		              sens(k,run) = -r_out'*X*r_in;  % compute sensitivity 
		              log_sens(k,run) = abs(xi*(-1)*r_out'*X*r_in)/(err(k,1)); % compute log-sensitivity           
	                end
           
		           log_sens(k,num+1) = mean(log_sens(k,1:num)); % compute norm of log-sensitivity 
		           sens(k,num+1) = mean(sens(k,1:num));
		       end
		       % sort results by nominal error and check consistency between err and err_calc     

		       Z = [err err_calc log_sens];
		       Z = sortrows(Z);
		       err = Z(:,1);
		       Z(:,1) = [];
		       err_calc = Z(:,1);
		       Z(:,1) = [];
		       log_sens = Z(:,1:num+1);
		       Diff = err-err_calc

		       % save results 
		       savetag = sprintf('../results/log_sens-rings/log_sens_%s_%d-ring_1-%d',opt,N,out);
		       save(savetag,'log_sens','sens','err','err_calc','Diff');
			end
	    end
	end
end

clear all; close all;
