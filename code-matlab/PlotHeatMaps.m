%% PlotHeatMaps.m
function PlotHeatMaps(density,Ctrl)

% SPDX-FileCopyrightText: Copyright (C) 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2022 S Shermer <lw1660@gmail.com>, Swansea University
% SPDX-FileCopyrightText: Copyright (C) 2023 Sean P O'Neil <seanonei@usc.edu>, University of Southern California
% 
% SPDX-License-Identifier: CC-BY-SA-4.0

% PlotHeatMaps(density,Ctrl) plots error density for chains 

% Input
%   density is a structure generated by 
%   density = AnalyzeErrorDensity(Err,Ctrl)
%   Ctrl - controller data

% Output
%   plot - plots depicting the error density

clf; set(gcf,'Position',[1600,800,1000,600])

%% Colormap

%% Histogram equalisation like colormap to maximise distinction
%cm = inferno(2^14);
%m = size(cm,1);
%buckets = histcounts(reshape(Density.^(1/10),numel(Density),1),m);
%buckets = floor(buckets / sum(buckets(2:end)) * (m-1));
%cm_eq(1,:) = cm(1,:); % Null bucket color
%idx = 2;
%for c = 2:m
%  cm_eq(c,:) = cm(idx,:);
%  idx = idx + buckets(c);
%end
cm_eq = inferno(2^14);

%% Plot Density

s_n = size(density.delta_kde,1);
bins = numel(density.delta_xi);
% note that median seems needed due to numerical rounding errors with diff
bin_size = median(diff(density.delta_xi));

subplot(3,3,[1:2, 4:5]);
DenseImg = (density.delta_kde.^(1/10))';
DenseImg = DenseImg(size(DenseImg,1):-1:1,:);
haxes = gca;
imagesc(haxes, DenseImg);
c = colorbar();
colourRange = caxis(haxes);
c.Ticks = linspace(colourRange(1),colourRange(end), 5);
c.TickLabels = arrayfun(@(n)sprintf('%0.2g', c.Ticks(n).^10),[1:5],'UniformOutput',false);
% c.Ticks.^10;

c.Label.String = 'Density';
title('Error density estimation');
xlabel('Strength');
s_steps = max(floor( (s_n-1) / 10),1);
xticks([0:s_steps:s_n-1]);
xticklabels([0:s_steps:s_n-1]/(s_n-1));
ylabel('Error');
steps = ceil(bins/10);
yticks([0:steps:bins]);
yticklabels(round(100*(density.max_err:-bin_size*steps:density.min_err))/100);
colormap(gca,cm_eq);

%% Plot Mean and Std. Deviation

hold on;
MeanImg = 1+(density.max_err-density.mean)*bins/(density.max_err-density.min_err);
StdDevImg = sqrt(density.var)*bins/(density.max_err-density.min_err);
plot(1:s_n,MeanImg,'Linewidth',1,'Color',[0 1 0]);
plot(1:s_n,MeanImg + StdDevImg,'-','Linewidth',1,'Color',[0 1 1]);
plot(1:s_n,MeanImg - StdDevImg,'-','Linewidth',1,'Color',[0 1 1]);

l = legend('Mean Error','Std. Dev.','Location','NorthWest');
l.Color = 'none';
l.TextColor = 'yellow';
l.EdgeColor = 'red';

title(sprintf('Error Density; Err(0) = %0.4g, log-sensitivity: %0.4g',density.mean(1),density.sensitivity));

%% Plot total distributions

subplot(3,3,[3,6]);

plot(density.total_kde,density.total_xi);
hold on; set(gca,'YAxisLocation','right');
line([0 max(density.total_kde)],density.total_mean*[1 1],'color','red');
hold off;

title(sprintf('Total Distribution\nMean: %0.4g, Variance: %0.4g', density.total_mean, density.total_var));
ylabel('Error');
xlabel('Density');
axis tight; 

%% Controller

subplot(3,3,7);

bar(Ctrl.x);
if isfield(Ctrl,'overlap')
  title(sprintf('T =%.4g; Overlap=%.4g',Ctrl.T,Ctrl.overlap));
else
  title(sprintf('T =%.4g',Ctrl.T));
end
xlabel('Spin');
ylabel('Bias');

subplot(3,3,8);
heatmaptext(round(100*Ctrl.rho_In_V)/100,'PRECISION',2,'COLORBAR',false,'CMAP',parula);
title('\rho_{IN}');

subplot(3,3,9);
heatmaptext(round(100*Ctrl.rho_Out_V)/100,'PRECISION',2,'COLORBAR',false,'CMAP',parula);
title('\rho_{OUT}');
%}
