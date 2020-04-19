%% Bias corrections factors for air temperature in GC
% Bryce Mihalevich
% 4/12/20

% Following methods from Teutschbein et al., 2012

%% Bias correction
% import era5 data 
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\wsx_comparison\era5Land_wsx.mat');

% go to directory containing figures comparing era5 data to sites in GC
cd('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\wsx_comparison\Figures_noPage')
k = dir('*AirTemp.fig');
% hand pick weather stations to use (GCMRC sites)
gcmrc_sites = [2 3 8 12 13 17 18 19 20 24]; %(3 includes page, consider not including or removing page from the analysis)

wsxGridID = fieldnames(era5_wsx.selected);
% loop through figures
for i = 1:length(gcmrc_sites)
    %% get data from the figure
    % find figure
    figInd = find(contains({k.name}',wsxGridID(gcmrc_sites(i)))); %index location of figure in k;
    % open figure
    uiopen(fullfile(k(figInd).folder,k(figInd).name),1);    
    
    af = gcf;
    h = af.Children(2).Children;
    dt = h(5).XData; %datetimes
    obsData = h(5).YData; % observed
    modData = h(4).YData; % era5
    % h(4).Color = [.5 .5 .5];
    close all %close the figure
    
    %% get monthly linear offset values
    % monthly average temperature for observed and gcm
    monthAvg_obsData = accumarray(dt.Month',obsData,[],@mean);
    monthAvg_modData = accumarray(dt.Month',modData,[],@mean);
    linearOffset(:,i) = monthAvg_obsData-monthAvg_modData;
    linearCorrection = modData' + linearOffset(dt.Month,i);
    
    %% get monthly variance multiplyer values
    % monthly standard deviation for observed and gcm
    linearCorrectionMonthly = accumarray(dt.Month',linearCorrection,[],@mean);
    zeroMeanSeries = linearCorrection - linearCorrectionMonthly(dt.Month);
    monthSTD_obsData = accumarray(dt.Month',obsData,[],@std);
    monthSTD_zeroMeanSeries = accumarray(dt.Month',zeroMeanSeries,[],@std);    
    varianceMultiplyer(:,i) = monthSTD_obsData./monthSTD_zeroMeanSeries;
end
%% save factors 
cd('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\bias_correction')
save('era5_airTempBiasCorrectionFactors_noPage.mat','linearOffset','varianceMultiplyer')
