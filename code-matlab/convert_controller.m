%% Convert ring controller data from .csv files to .mat (convert_controller.m)

% SPDX-FileCopyrightText: Copyright (C) 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2022 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2023 S Sean Patrick O'Neil
% <seanonei@usc.edu>, University of Southern California 
% 
% SPDX-License-Identifier: CC-BY-SA-4.0

function sys = convert_controller(data,in,out,optimiser,N,x)

  % Converts controller data in .csv format to usable .mat format 
  
  % This takes the following array in input (requires use of readmatrix()
  % to convert the .csv files in the ../data directory to an array:
  % array columns

  %   1 to N - bias values 
  %   N+2    - fidelity
  %   N+1    - Time

  % Outputs a cell array sys with each entry a structure s. that contains
  % the fidelity error, time of transfer, and bias fields
  
  switch x 
      case 1
        H = zeros(N,N);
        H = diag(ones(1,N-1),1);
        H = H+H';
      case 2
        H = zeros(N,N);
        H = diag(ones(1,N-1),1);
        H = H+diag(ones(1,1),N-1);
        H = H+H';
  end
 
  rho0 = zeros(N);  rho0(in,in)   = 1;
  rho1 = zeros(N);  rho1(out,out) = 1;

  sys = cell(1,size(data,1));
  for l = 1:size(data,1)
    s.Nspin = N;
    s.In = in;
    s.Out = out;
    s.Optimiser = optimiser;
    s.fidelity = data(l,N+2);
    s.T = data(l,N+1);
    s.x = data(l,1:N);
    s.H = H + diag(s.x);

    % Calculate steady state info
    [V,e] = eig(s.H);
    s.H_V = V;
    s.H_e = diag(e)';
    s.rho_In_V = V'*rho0*V;
    s.rho_Out_V = V'*rho1*V;
    s.steady_state = diag(s.rho_In_V)';
    s.purity = sum(s.steady_state.^2);

    sys{l} = s;
  end

end