%% Produce basis matrices (Gell-Mann matrices) for Bloch formalism (bloch_basis.m)

function [b]=bloch_basis(N,I)

% SPDX-FileCopyrightText: Copyright (C) 2023 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-FileCopyrightText: Copyright (C) 2023 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2023 SM Shermer <lw1660@gmail.com>
% SPDX-License-Identifier: CC-BY-SA-4.0 

% This function takes as input the size of the system (N) and generates the
% N^2 Gell-Mann matrices a a basis for the Hermitian matrices 

% Produce Natural Basis Vectors
for j=1:N
    e{j}=I(:,j);
end

% Produce Diagonal Basis Matrices 
for k=1:(N-1)
    b{k}=zeros(N,N);
    for l=1:k
        b{k}=b{k}+e{l}*e{l}';
    end
    b{k}=1/sqrt(k+k^2)*(b{k}-k*e{k+1}*e{k+1}');
end

% Re-Index Diagonal Basis Matrices
for k=1:(N-1)
    b{N^2-N+k}=b{k};
end

b{N^2}=(1/sqrt(N))*eye(N,N);

% Produce Real Symmetric Off-Diagonal Basis Matrices 
z=1;
for k=1:N
    for l=(k):(N-1)
        b{z}=1/sqrt(2)*(e{k}*e{l+1}'+e{l+1}*e{k}');
        z=z+1;
    end
end

% Produce Imaginary Hermitian Off-Diagonal Basis Matrices
for k=1:N
    for l=(k):(N-1)
        b{z}=1/sqrt(2)*(i*e{k}*e{l+1}'-i*e{l+1}*e{k}');
        z=z+1;
    end
end

% Test Orthogonality of Basis Matrices
for j=1:N^2
    for k=1:N^2
        O(j,k)=trace(b{j}*b{k});
    end
end

