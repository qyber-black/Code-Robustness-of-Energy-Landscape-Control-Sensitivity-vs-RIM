%% Calcualte error data for kde robusntess calculations (AnalyzeRobustnessDephasing.m)

% SPDX-FileCopyrightText: Copyright (C) 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2022 S Shermer <lw1660@gmail.com>, Swansea University
% 
% SPDX-License-Identifier: CC-BY-SA-4.0

function [err,r] = AnalyzeRobustnessDephasing(sys,G,in,out)
%
% function [err,r] = AnalyzeRobustness(sys,G,in,out) analyses
% the robustness of the controllers for system sys with regard to
% the dephasing operators G for transfer from state in to out
%
% input:  sys -- structure with fields
%         H   -- total Hamiltonian (including biases)
%         T   -- transfer time T
%         G   -- list of dephasing operators
%         in  -- index of input node
%         out -- index of output node
%
% output: err -- vector of distances of output r{k} from target r1
%         r   -- list of output vectors for each dephasing op
%

H = sys.H;
T = sys.T;
if isfield(sys,'dT') % if windowed readout
  windowed = 1;
  dT = sys.dT;
else
  windowed = 0;
end

[v,e] = eig(H);
e = diag(e);
N = length(e);

r0 = zeros(N); r0(in,in)  =1; r0 = v'*r0*v;
r1 = zeros(N); r1(out,out)=1; r1 = v'*r1*v;

omega = e*ones(1,N)-ones(N,1)*(e');
Gamma = G+G';
Gamma = Gamma/sum(Gamma(1:end));
gamma = [0:1000]*(0.1/1000);

for k = 1:length(gamma)
  x = -1i*omega - gamma(k)*Gamma;
  r{k} = r0.*exp(T*x);
  err(k)  = 1 - real(trace(r{k}'*r1));
end