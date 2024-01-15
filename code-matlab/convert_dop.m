%% Convert dephasing operator data from .csv to .mat format (convert_dop.m)

% SPDX-FileCopyrightText: Copyright (C) 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2022 S Shermer <lw1660@gmail.com>, Swansea University
% 
% SPDX-License-Identifier: CC-BY-SA-4.0

function dop = convert_dop(arr,N,test)
% Creates triangular matrices for py dephasing ops (optionally test and 
% plot 2D projections)
%
% Input:
%   arr - array of dephasing operators taken from .csv file in the ../data
%         directory
%   N - size of ring 
%   test - option to test physical realization of dephasing operators (0 =
%          no test / 1 = test)
%   plot - option to plot dephasing operators (0 = no plot / 1 = plot) 
%
% Outupt:
%   dop - 100,000 element cell array, with each cell containing an NxN
%         lower-triangular matrix of (normalized) dephasing operators  

  dop = cell(1,size(arr,1));
  for l = 1:size(arr,1)
    g = triu(ones(N),1);
    g(g~=0) = arr(l,:);
    g = g';

    if (test > 0)
      [VG,state] = InvertDephasing(g);
    end
    
    dop{l} = g;
end