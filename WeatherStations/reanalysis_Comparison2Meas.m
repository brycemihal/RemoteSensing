%% Script to compare reanalysis data to measured data
% Bryce Mihalevich
% 1/24/2020

% info
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview

% Relative path to files
n = mfilename;
p = mfilename('fullpath');
cd(p(1:end-length(n)))
mfileDir = pwd;

%% Get data from database
run('queryWeatherData_LowerBasinWSX.m')

% fix some bad data
loc = find(cell2mat(lowerBasinWSX.Mile_24_5_Upper.Data.AirTemp(:,2))<-20);
lowerBasinWSX.Mile_24_5_Upper.Data.AirTemp(loc,:) = [];

% remove page
lowerBasinWSX = rmfield(lowerBasinWSX,'Page_Municipal_Airport');
%% Import data from the netcdf
ncFileDir = 'C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\raw_data';
if ~exist('era5Land_wsx_noPage.mat','file')>0
    cd(ncFileDir)
    k = dir('*.nc');
    count = 0;
    % loop for each file. Files should be ordered by ascending years (05,06,..etc)
    for i = 1:length(k)
        netCDF_file = k(i).name;
        %     ncinfo(netCDF_file);
        era5.air_T = squeeze(ncread(netCDF_file,'t2m'));
        era5.dewPoint_T = squeeze(ncread(netCDF_file,'d2m'));
        era5.solar_down = squeeze(ncread(netCDF_file,'ssrd')); %this variable is the model equivalent of what would be measured by a pyranometer
        era5.wind_U = squeeze(ncread(netCDF_file,'u10'));
        era5.wind_V = squeeze(ncread(netCDF_file,'v10'));
        era5.latitude = (ncread(netCDF_file,'latitude'));
        era5.longitude = (ncread(netCDF_file,'longitude'));
        era5.time_utc = datetime('1900-01-01 00:00:00') + hours(ncread(netCDF_file,'time'));
        era5.time_mst = era5.time_utc-hours(7);
        
        %% Import shapefile data
        USshape = shaperead('C:\Users\Echo\Documents\GIS\cb_2016_us_state_5m\cb_2016_us_state_5m.shp');
        %     lakePowell = shaperead('C:\Users\Echo\Documents\ColoradoRiverBasin\Weather Data\Reanalysis\ERA5Land\code\shp\LakePowell_wgs.shp');
        %     lakeMead = shaperead('C:\Users\Echo\Documents\ColoradoRiverBasin\Weather Data\Reanalysis\ERA5Land\code\shp\LakeMead_wgs.shp');
        gc_centerline = shaperead('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\grids_to_model_elements\shp\GrandCanyon_CenterLine_wgs.shp');
        gc_centerline_abvLF = shaperead('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\grids_to_model_elements\shp\GrandCanyon_CenterLine_abvLF_wgs.shp');
        
        load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\grids_to_model_elements\lakeMeadPolygon.mat')
        load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\grids_to_model_elements\lakePowellPolygon.mat')
        
        % Import river element junction data
        riverFile = 'ColoradoRiver_GrandCanyon_junctions.csv';
        fid = fopen(riverFile,'r');
        riverData = textscan(fid,'%s%f%f%f%f%[^\n\r]','HeaderLines',1,'Delimiter',',');
        fclose(fid);
        
        rivNodeX = riverData{2};
        rivNodeY = riverData{3};
        rivNodeKM = riverData{5};
        
        for j = 1:length(riverData{1})-1
            % calculate centriod xy
            rivCentroidX(j) = mean(rivNodeX(j:j+1));
            rivCentroidY(j) = mean(rivNodeY(j:j+1));
            
            % determine centriod river km
            rivCentroidKM = rivNodeKM(1:end-1) + round(diff(rivNodeKM),1)/2;
        end
        % convert UTM to lat lon
        [rivCentroid_lat, rivCentroid_lon] = utmups_inv(rivCentroidX,rivCentroidY, 12, 1);
        
        % lakePowell = shaperead('C:\Users\Echo\Documents\ColoradoRiverBasin\Weather Data\Reanalysis\ERA5Land\code\shp\LakePowell_wgs.shp');
        % lakeMead = shaperead('C:\Users\Echo\Documents\ColoradoRiverBasin\Weather Data\Reanalysis\ERA5Land\code\shp\LakeMead_wgs.shp');
        
        %     lp_pgon = polyshape(lakePowell(1).X,lakePowell(1).Y);
        %     lm_pgon = polyshape(lakeMead.X,lakeMead.Y);
        %% Determine what reanalysis grids correspond to weather stations
        %get all unique pairs of grid points
        [era5_lat,era5_lon] = meshgrid(era5_wsx.latitude,era5_wsx.longitude);
        era5_lat = era5_lat(:);
        era5_lon = era5_lon(:);
        
        sites = fieldnames(lowerBasinWSX);
        for j = 1:length(sites)
            % caclulate distance between weather station and RS cells and grab index of minimum distances
            site_lat(j) = lowerBasinWSX.(sites{j}).Latitude;
            site_lon(j) = lowerBasinWSX.(sites{j}).Longitude;
            [~,ind(j)] = min(sqrt((lowerBasinWSX.(sites{j}).Longitude-era5_lon).^2 + (lowerBasinWSX.(sites{j}).Latitude-era5_lat).^2));
        end
        
        [uInd,~,siteSubs] = unique(ind);
        
        %%
        close all
        figure()
        hold on        
        % add era5 grid boxes
        lat_spacing = mode(diff(unique(era5_lat)))/2;
        lon_spacing = mode(diff(unique(era5_lon)))/2;
        lat_lines = [min(era5_lat)-lat_spacing; unique(era5_lat)+lat_spacing];
        lon_lines = [min(era5_lon)-lon_spacing; unique(era5_lon)+lon_spacing];
        for j = 1:length(lat_lines)
            plot([min(lon_lines) max(lon_lines)],[lat_lines(j) lat_lines(j)],'k')
        end
        for j = 1:length(lon_lines)
            plot([lon_lines(j),lon_lines(j)],[min(lat_lines) max(lat_lines)],'k')
        end
        
        % add patch for selected boxes
        for j = 1:length(uInd)
            xq = [era5_lon(uInd(j))-lon_spacing era5_lon(uInd(j))+lon_spacing...
                era5_lon(uInd(j))+lon_spacing era5_lon(uInd(j))-lon_spacing]; 
            yq = [era5_lat(uInd(j))-lat_spacing era5_lat(uInd(j))-lat_spacing...
                era5_lat(uInd(j))+lat_spacing era5_lat(uInd(j))+lat_spacing];
            p = patch(xq,yq,'red');
            p.FaceColor = [1 0 0];
            p.FaceAlpha = .5;
            t = text();
            t.String = num2str(j);
            t.Units = 'data';
            t.FontSize = 10;
            t.FontWeight = 'bold';
            t.Position = [era5_lon(uInd(j)) era5_lat(uInd(j)) 0];
            t.HorizontalAlignment = 'center';
            t.VerticalAlignment = 'middle';
            t.Clipping = 'on';
        end

        % add patch for gcmrc boxes
        gcmrc_sites = [2 3 8 12 13 17 18 19 20 24];
        for j = 1:length(gcmrc_sites)
            xq = [era5_lon(uInd(gcmrc_sites(j)))-lon_spacing era5_lon(uInd(gcmrc_sites(j)))+lon_spacing...
                 era5_lon(uInd(gcmrc_sites(j)))+lon_spacing era5_lon(uInd(gcmrc_sites(j)))-lon_spacing]; 
            yq = [era5_lat(uInd(gcmrc_sites(j)))-lat_spacing era5_lat(uInd(gcmrc_sites(j)))-lat_spacing...
                era5_lat(uInd(gcmrc_sites(j)))+lat_spacing era5_lat(uInd(gcmrc_sites(j)))+lat_spacing];
            p = patch(xq,yq,'red');
            p.FaceColor = [1 0 0];
            p.FaceAlpha = .5;
        end                
        
        plot(rivCentroid_lon,rivCentroid_lat,'-','LineWidth', 2, 'Color',[0.5961 0.8667 0.9882])
        gcmrc_sites2 = [15, 20:33];
        plot(site_lon(gcmrc_sites2),site_lat(gcmrc_sites2),'ok','MarkerFaceColor','k')
        plot_google_map('MapType','terrain','ShowLabels',0)
        ax = gca;
        ax.YLabel.String = 'Latitude';
        ax.XLabel.String = 'Longitude';
        ax.XLim = [-113.7609 -111.4235];
        ax.YLim = [35.6008   37.0857];
%         ax.TickDir = 'out';
%         ax.XMinorTick ="on";
%         ax.YMinorTick ="on";
%         ax.Box = 'on';
%         
        str = ['grid = ERA dataset' newline 'black dots = weather stations' newline 'red boxes = nearest ERA grid'];
        t = text();
        t.String = str;
        t.Units = 'normalized';
        t.Position = [0.01 .9 0];
        t.FontSize = 12;
        t.BackgroundColor = 'w';
        
        %% only plot gcmrc grids
        % gcmrc_sites = [2 3 8 12 13 17 18 19 20 24];

        %% Get data from the idendified grid indexes
        for j = 1:length(uInd)
            
            % get era5 lat and lon index
            era5_lat_ind = find(ismember(era5.latitude,era5_lat(uInd(j))));
            era5_lon_ind = find(ismember(era5.longitude,era5_lon(uInd(j))));
            id = strcat('r',num2str(era5_lon_ind),'c',num2str(era5_lat_ind));
            
            %% air temperature
            qAirT = squeeze(cell2mat({era5.air_T(era5_lon_ind,era5_lat_ind,:)-273.15}));
            
            if i == 1
                era5.selected.(id).AirTemp = {era5.time_mst qAirT};
            else % append additional years
                era5.selected.(id).AirTemp{1} = [era5.selected.(id).AirTemp{1}; era5.time_mst];
                era5.selected.(id).AirTemp{2} = [era5.selected.(id).AirTemp{2}; qAirT];
            end
            clear qAirT
            
            %% wind
            wU = squeeze(cell2mat({era5.wind_U(era5_lon_ind,era5_lat_ind,:)}));
            wV = squeeze(cell2mat({era5.wind_V(era5_lon_ind,era5_lat_ind,:)}));
            qWind = sqrt(wU.^2+wV.^2);
            
            if i == 1
                era5.selected.(id).WindSpeed = {era5.time_mst qWind};
            else % append additional years
                era5.selected.(id).WindSpeed{1} = [era5.selected.(id).WindSpeed{1}; era5.time_mst];
                era5.selected.(id).WindSpeed{2} = [era5.selected.(id).WindSpeed{2}; qWind];
            end
            clear qWind
            
            %% solar
            qSSRD = squeeze(cell2mat({era5.solar_down(era5_lon_ind,era5_lat_ind,:)}))/3600;
            
            % loop through each timestep
            qSolar(1,1) = 0;
            for jj = 2:length(era5.time_utc)
                if hour(era5.time_utc(jj))==1
                    qSolar(jj,1) = qSSRD(jj);
                else
                    qSolar(jj,1) = qSSRD(jj)- qSSRD(jj-1);
                end
            end
            
            if i == 1
                era5.selected.(id).SolarRadiation = {era5.time_mst qSolar};
            else % append additional years
                era5.selected.(id).SolarRadiation{1} = [era5.selected.(id).SolarRadiation{1}; era5.time_mst];
                era5.selected.(id).SolarRadiation{2} = [era5.selected.(id).SolarRadiation{2}; qSolar];
            end
            clear qSolar
            
            %% relative humidity
            qDewPT = squeeze(cell2mat({era5.dewPoint_T(era5_lon_ind,era5_lat_ind,:)-273.15}));
            qAirT = squeeze(cell2mat({era5.air_T(era5_lon_ind,era5_lat_ind,:)-273.15}));
            esat = 4.596 * exp(17.27*qAirT./(237.3+qAirT));
            eair = 4.596 * exp(17.27*qDewPT./(237.3+qDewPT));
            qRH = 100*eair./esat;
            
            if i == 1
                era5.selected.(id).RelativeHumidity = {era5.time_mst qRH};
            else % append additional years
                era5.selected.(id).RelativeHumidity{1} = [era5.selected.(id).RelativeHumidity{1}; era5.time_mst];
                era5.selected.(id).RelativeHumidity{2} = [era5.selected.(id).RelativeHumidity{2}; qRH];
            end
            clear qRH
            
            count = count+1;
            if mod(count,round((length(k)*length(uInd))/100))==0 %write out percent progress to command window
                percentComplete = round(count/(length(k)*length(uInd))*100);
                pp = sprintf('~%d%% complete',percentComplete);
                disp(pp)
            end
        end
    end
    era5_wsx = era5;
    uInd_wsx = uInd;
    subs_wsx = siteSubs;
    cd(mfileDir)
    save('era5Land_wsx_noPage.mat','era5_wsx','uInd_wsx','subs_wsx')
else
    load('era5Land_wsx.mat')
end
%% Calculate statistics comparing weather station to RS grid
% loop for each grid
params = {'AirTemp','RelativeHumidity','WindSpeed','SolarRadiation'};
units = {'[C]','[%]','[m/s]','[W/m2]'};
gridID = fieldnames(era5.selected);

for g = 1:length(uInd)
    % determine corresponding weather stations for grid
    loc = find(siteSubs==g); %index location of sites within grid
    
    for i =1:length(sites(loc))
        gridSites{g,i} = char(sites(loc(i)));
    end
    
    for p = 1:length(params)
        % assign modeled data
        mod_times = era5.selected.(gridID{g}).(params{p}){1};
        mod_values = era5.selected.(gridID{g}).(params{p}){2};
    
        [E{g,p},RMSE{g,p},monthlyRMSE{g,p},minE{g,p},maxE{g,p},ME{g,p},MAE{g,p}]...
            = calculateStats(mod_times,mod_values,lowerBasinWSX,sites,params{p},units{p},loc,gridID{g});
    end
end

save('reanalysisStats.mat','E','RMSE','monthlyRMSE','minE','maxE','ME','MAE','gridID','params')


%% monthly rmse
for i = 1:size(monthlyRMSE,1)
    for j = 1:size(monthlyRMSE,2)
        if isnan(monthlyRMSE{i,j})
            monthlyRMSE{i,j} = NaN(12,1);
        end
    end
end
monthlyRMSE{6,2}(1:12) = NaN(12,1);
%% Plot monthly RMSE
monthStr = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for i = 1:size(monthlyRMSE,2)
    pData = [monthlyRMSE{:,i}];
    pData(pData<=0) = NaN;
    pDataAvg = nanmean(pData,2);
    
    figure()
    hold on
    for j = 1:size(pData,2)
        h1 = plot(pData(:,j),'Color',[.5 .5 .5]);
    end
    h2 = plot(pDataAvg,'LineWidth',2,'Color','k');
    
    ax = gca;
    ax.XLim = [1 12];
    ax.XTickLabel = monthStr;
    ax.YLabel.String = params{i};
    l = legend([h1 h2],{'All grids','Average'});
end

%% Linear model
uiopen('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\wsx_comparison\Figures\ERA5r23c9_AirTemp.fig',1);
af = gcf;
h = af.Children(2).Children;
obsData = h(5).YData; % observed
modData = h(4).YData; % era5

l_mdl = fitlm(obsData,modData);
figure()
plot(l_mdl)

%% function to calculate statistics
function[E,RMSE,monthlyRMSE,minE,maxE,ME,MAE] = calculateStats(mod_times,mod_values,lowerBasinWSX,sites,param,units,loc,gridID)

% combine observed data among sites that fall within grid
combined_times = [];
combined_values = [];
for s = loc'
    try
        obs_times = datetime(lowerBasinWSX.(sites{s}).Data.(param)(:,1)); %observed data
    catch
        E = NaN;
        RMSE = NaN;
        minE = NaN;
        maxE = NaN;
        ME = NaN;
        MAE = NaN;
        monthlyRMSE =  NaN;
        return
    end
    obs_values = cell2mat(lowerBasinWSX.(sites{s}).Data.(param)(:,2)); %observed datetimes
    
    combined_times = [combined_times; obs_times]; %stack times
    combined_values = [combined_values; obs_values]; %stack values
end

[obs_times,ord] = sort(combined_times); %sort by datetime
obs_values = combined_values(ord); % apply sort order to values

% average observed data to houly interval
obs_times_str = datestr(obs_times,'mm/dd/yyyy HH'); % convert to string to only show the hour
[uDateHour,~,subs] = unique(obs_times_str,'stable','rows'); % get index of unique date-hours

obs_times_hourly = datetime(uDateHour,'InputFormat','MM/dd/yyyy HH'); %revert back to datetime format
obs_values_hourly = accumarray(subs,obs_values,[],@mean); % average within hour

% match observed data to era5 data
[common_times,obsInd,modInd] = intersect(obs_times_hourly,mod_times); % times common to both
Qobs = obs_values_hourly(obsInd); % values corresponding to common times
Qmod = mod_values(modInd);  % values corresponding to common times

% caclualte statistics and save them to cell {grid x parameter}
rmse_fun = @(E) sqrt(nanmean(E.^2)); %rmse function

E = Qobs-Qmod; %hourly error
RMSE = rmse_fun(E); %rmse
minE = min(abs(E)); %minimum abosolute error
maxE = max(abs(E)); % maximum abosolute error
ME = nanmean(E); % mean abosolute error
MAE = nanmean(abs(E)); % mean error
monthlyRMSE = accumarray(common_times.Month,E,[],rmse_fun);

%send to function for plotting;
plotData(common_times,Qobs,Qmod,gridID,sites(loc),param,units)
end

%% Save figures to word document
function[] = savePlots2Word()
save2Dir = 'C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\wsx_comparison\';
cd(save2Dir)
rptName = 'ERA5_reanalysis_lowerBasinWSX_comparison_plots_noPage.doc';
generateReport(rptName)
end

%% Save all figures as the Matlab .fig format so that they can be reopend
function[] = savePlotsAsFig(site,parameter)
save2Dir = 'C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\wsx_comparison\';
cd(save2Dir)
figDir = fullfile(save2Dir,'Figures_noPage');
if ~exist(figDir,'dir')
    mkdir(figDir)
end
cd(figDir)

rptName = strcat(site,'_',parameter);
figs = findobj('Type','Figure');
for i = 1:length(figs)
    figName = rptName;
    savefig(figs(i),figName);
end
end

%% Plot Function
function[] = plotData(common_times,Qobs,Qmod,gridID,sites,params,units)
rmse_fun = @(E) sqrt(nanmean(E.^2)); %rmse function
E = Qobs-Qmod;
RMSE = rmse_fun(E);

figure()
set(gcf,'Name',strcat('ERA5',gridID,'_',char(params)));
subplot(2,2,4)
hold on
h = histogram(E);
h.Orientation = 'horizontal';
h.FaceColor = [.5 .5 .5];
h.EdgeColor = 'none';
yline(0,'Color','k');
ax3 = gca;
ax3.YTickLabel = [];
ax3.Position = [.83 .08 .16 .25];

subplot(2,2,3)
hold on
plot(common_times,E,'Color',[.5 .5 .5])
yline(0,'Color','k');
ax2 = gca;
ax2.Position = [0.1 .08 .72 .25];
ax3.YLim = ax2.YLim;
ax2.YLabel.String = 'Residual';

subplot(2,2,1:2)
hold on
plot(common_times,Qobs,'k')
plot(common_times,Qmod,'r--')
ax1 = gca;
ax1.XLim(1) = common_times(1);
ax1.XLim(2) = common_times(end);
ax1.Position = [0.1 .36 .72 .60];
ax1.YLabel.String = strcat((params),{' '},units);
ax1.XTickLabel = [];
l = legend('Observed','ERA5');
l.Orientation = 'horizontal';
l.Position(1) = .11;
l.Position(2) = .93;
l.Position(4) = .05;

t = text();
t.String = sprintf(strcat('RMSE:\n',num2str(RMSE)));
t.Units = 'normalized';
t.Position = [1.02 .1 0];

t2 = text();
t2.Interpreter = 'none';
str = 'Sites:\n';
for i = 1:length(sites)
    str = strcat(str,sites{i},'\n');
end
t2.String = sprintf(str);
t2.Units = 'normalized';
t2.Position = [1.02 .9 0];
t2.FontSize = 8;

t3 = text();
t3.Units = 'normalized';
t3.String = sprintf(strcat('ERA5 Grid ID:\n',gridID));
t3.Position = [1.02 .5 0];
t3.FontSize = 9;

ax2.XLim = ax1.XLim;
savePlotsAsFig(strcat('ERA5',gridID),char(params))
savePlots2Word()
close all;
end
