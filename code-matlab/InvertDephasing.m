%% Test dephasing operators for physical realizability (InvertDephasing.m)

% SPDX-FileCopyrightText: Copyright (C) 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2022 S Shermer <lw1660@gmail.com>, Swansea University
% 
% SPDX-License-Identifier: CC-BY-SA-4.0

function [V,state,c,Gamma2] = InvertDephasing(G,domg)
% Invert dephasingRate and test for physical realization 
%
% Input:
%   G - lower-triangular matrix of dephasing operators 
%
% Output:
%   V - matrix of "inverted" dephasing operators used to test 
%   state - level of constraint violation if violation occurs 

tol = 1e-12;
state = 0;
N = size(G,1);
if ~exist('domg','var');
    domg = zeros(N);
end
m0 = 2;
V = zeros(N,N-1);

V(2,1) = sqrt(2*G(2,1));
if V(2,1)<tol
    m0 = m0+1;
end

for n=3:N
    for m=1:n-m0
        V(n,m) = ( G(n,1)+G(m+1,1)-G(n,m+1)+1i*domg(m+1,n) - sum(V(n,:).*conj(V(m+1,:))) ) / V(m+1,m);
    end
    const = 2*G(n,1)-sum(abs(V(n,:)).^2);
    if const>tol
        V(n,n-m0+1) = sqrt(const);
    elseif const<-tol
        disp(sprintf('Constraint violation at level %i',n))
        state=n; return
    else
      m0 = m0+1;
    end
end


end
