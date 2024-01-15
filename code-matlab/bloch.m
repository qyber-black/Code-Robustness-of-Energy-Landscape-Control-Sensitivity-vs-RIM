%% Produce Bloch-transformed matrices (bloch.m)

function [A,lambda,V,r_in,r_out] = bloch(N,Hd,r0,r1,b)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function takes as input the controlled Hamiltonian Hd, the input and output 
% states, and the Gell-Mann matrices and on output provides the
% real representation of Hd (A) in the Bloch formalism along with the
% vectors r_in and r_out

% Define input and output 
[V,lambda] = eig(Hd);
r1 = V'*r1*V;
r0 = V'*r0*V;

% Produce Bloch Representations of Hd (A)
for n=1:N^2
    for m=1:N^2
        A(n,m)=trace(1*i*lambda*(b{n}*b{m}-b{m}*b{n}));
    end
    r_in(n,1) = trace(b{n}*r0);
    r_out(n,1) = trace(b{n}*r1);
end


