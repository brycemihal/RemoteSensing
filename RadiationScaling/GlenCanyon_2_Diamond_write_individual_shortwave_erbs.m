% Script to calculate diffuse radiation using Erbs et al method
clear; clc;

% cd to location of m file
m = mfilename;
p = mfilename('fullpath');
cd(p(1:end-length(m)))
mfileDir = pwd;

%% Import pre-calculated data
% ERA5
% includes:
% uInd: index values of grid cells that overlap river elements
% subs: index of what river kilometers correspond to specific grid cells
% era5: structure containing "selected" grid cells
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\grids_to_model_elements\era5Land.mat')

% get all unique pairs of grid points
[era5_lat,era5_lon] = meshgrid(era5.latitude,era5.longitude);
era5_lat = era5_lat(:);
era5_lon = era5_lon(:);

% Solar incidence with Leap Day
% includes:
% reachAvgIncidence: 1 year of incidence data between GCD and spencer creek at 15 minute increments
% dateTimes: datetimes array corresponding to incidence values
% incidence: incidence every 100m (not used)
% skyVF: SVF every 100m (not used)
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\radiation_scaling\SolarIncidence_output_366.mat') % for datetimes
load('reachAvgIncidence_GCD_to_SpencerCk_wLeapDay.mat')
dateTimes_LY = dateTimes; % march 1 repeats twice
reachAvgIncidence_LY = reachAvgIncidence;

% Solar incidence data
% includes:
% reachAvgIncidence: 1 year of incidence data between GCD and spencer creek at 15 minute increments
% dateTimes: datetimes array corresponding to incidence values (Does not include leap year)
% incidence: incidence every 100m (not used)
% skyVF: SVF every 100m (not used)
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\radiation_scaling\SolarIncidence_output_365.mat') % for datetimes
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\radiation_scaling\reachAvgIncidence_GCD_to_SpencerCk.mat') % for incidence
clear('incidence','skyVF')

% Sky view factor data
% includes:
% reachAvgSVF: SVF values averaged by model element between GCD and spencer creek
load('reachAvgSVF_GCD_to_SpencerCk.mat') % import svf da

% Model Link file (junction name and riv km)
filename = 'GrandCanyon_model_links_magirl.csv';
fid = fopen(filename,'r');
modLinks = textscan(fid,'%s%f%[^\n\r]','HeaderLines',1,'Delimiter',',');
fclose(fid);
elementRivKM = modLinks{2};
elementLinkNames = elementRivKM(1:end-1)+round((round(diff(elementRivKM),1)/2),2);

%% Remove elements downstream of Diamond Creek from incidence and svf data
writeElementLinkNames = elementLinkNames;
loc = find(writeElementLinkNames > 362.7,1);
reachAvgIncidence(:,loc+1:end) = [];
reachAvgIncidence_LY(:,loc+1:end) = [];
reachAvgSVF(loc+1:end) = [];
writeElementLinkNames(loc+1:end) = [];

%% reduce data down to hourly interval
dateTimesHourly = (dateTimes(1):hours(1):dateTimes(end))'; %hourly
reachAvgIncidence = interp1(dateTimes,reachAvgIncidence,dateTimesHourly);
dateTimes = dateTimesHourly;

dateTimes_LY = dateTimes_LY+calyears(1);
dateTimes_LY = (dateTimes_LY(1):minutes(15):dateTimes_LY(end))';
dateTimesHourly = (dateTimes_LY(1):hours(1):dateTimes_LY(end))'; %hourly
reachAvgIncidence_LY = interp1(dateTimes_LY,reachAvgIncidence_LY,dateTimesHourly);
dateTimes_LY = dateTimesHourly;

%% repeat for number of years in ERA5 dataset
minYear = year(min(era5.selected.r10c9.AirTemp{1}));
maxYear = year(max(era5.selected.r10c9.AirTemp{1}));

uYears = unique(year(dateTimes)); % dateTimes are 2019
if length(uYears)>1
    error('to many years in dateTimes array')
end

uYears_LY = unique(year(dateTimes_LY)); % dateTimes_LY are 2020
if length(uYears)>1
    error('to many years in dateTimes array')
end

% years to write out
writeYears = minYear:maxYear; %ascending order
writeDateTimes = [];
repReachAvgIncidence = [];
for i = 1:length(writeYears)
    if ~leapyear(writeYears(i))
        writeDateTimes = [writeDateTimes; dateTimes-calyears(abs(uYears - writeYears(i)))];
        repReachAvgIncidence = [repReachAvgIncidence; reachAvgIncidence];
    else
        writeDateTimes = [writeDateTimes; dateTimes_LY-calyears(abs(uYears_LY - writeYears(i)))];
        repReachAvgIncidence = [repReachAvgIncidence; reachAvgIncidence_LY];
    end
end

% repReachAvgIncidence = repmat(reachAvgIncidence,length(writeYears),1);

%% remove datetimes / incidence rows not in ERA5
eraFields = fieldnames(era5.selected);
era5_SolarTimes = era5.selected.(eraFields{1}).SolarRadiation{1};
ind = ismember(writeDateTimes,era5_SolarTimes);
writeDateTimesStr = datestr(writeDateTimes(ind),'mm/dd/yyyy HH:MM:SS');
repReachAvgIncidence = repReachAvgIncidence(ind,:);

%% Compute diffuse fraction
% loop through eraFields or uInd (same length)
matlabDataFile = 'era5_with_radiation_factors.mat';
if ~exist(matlabDataFile,'file')>0
    for i = 1:length(uInd)
        %% get solar data
        era5_Solar = era5.selected.(eraFields{i}).SolarRadiation{2}';
        era5_SolarTimes = era5.selected.(eraFields{i}).SolarRadiation{1};
        
        %% determine latitude and longitude of grid cell
        lat_deg = era5_lat(uInd(i));
        lon_deg = era5_lon(uInd(i));
        
        %% calculate top of atmosphere solar radiation
        DOYs = day(era5_SolarTimes,'dayofyear');
        LocalTimes = hour(era5_SolarTimes)+minute(era5_SolarTimes)/60;
        
        UTCoffset = -6.5; %timezone correction % should be -7 but for some reason -6.5 alings better.
        UTCtimes = LocalTimes - UTCoffset;
        
        for j = 1:length(DOYs)
            DOY = DOYs(j);
            [RsTOA(j),zenith_angle_deg(j),azimuth_angle_deg(j),sunrise,sunset,solar_decl,hour_angle]...
                = TOA_incoming_solar(DOY,UTCtimes(j),UTCoffset,lat_deg,lon_deg);
        end
        
        %% Check series
        %     close all
        %     hold on
        %     plot(commonTimes,RsTOA)
        %     plot(commonTimes,era5_Solar)
        %     pause()
        
        %% Calculate diffuse fractions following Erbs et al 1982 model
        kt = era5_Solar./(RsTOA.*sind(90-zenith_angle_deg));
        kt(kt>1)=1; %kt should not exceed 1.
        kt(isnan(kt))=0; %kt should not be nan
        kt = kt';
        % calculate diffuse radiation Above Canyon
        for j = 1:length(kt(:,1))
            if kt(j) > 0.8
                Jdiff_AC(j) = era5_Solar(j) * 0.165;
            elseif kt(j) >= 0.22 && kt(j) <= 0.8
                Jdiff_AC(j) = era5_Solar(j) * (0.9511-0.1604*kt(j)+4.38*kt(j)^2-16.64*kt(j)^3+12.34*kt(j)^4);
            elseif kt(j) >= 0 && kt(j) < 0.22
                Jdiff_AC(j) = era5_Solar(j) * (1-0.09*kt(j));
            end
        end
        
        %% apply diffuse fraction to subset up writeReachAvgIncidence
        gridRiverKMs = writeElementLinkNames(subs==i); % river km's that correspond to grid
        gridReachAvgSVF = reachAvgSVF(subs==i); % sky view factors that correspond to river km's
        gridReachAvgIncidence = repReachAvgIncidence(:,subs==i); % incidence factors that correspond to river km's
        
        %% calculate individual SW components and radiation factors for grid
        
        Jsn_interp = era5_Solar'; % Measured above canyon
        
        Jdiff_interp = Jdiff_AC'; % Diffuse above canyon
        
        Jdir_interp = Jsn_interp - Jdiff_interp; % Direct above canyon
        
        land_albedo = 0.01;
        
        Jreflected = (Jdiff_interp+Jdir_interp).*(1-gridReachAvgSVF')*land_albedo; %reflected in canyon
        Jdir =  Jdir_interp.*(1-gridReachAvgIncidence); %direct in canyon
        Jdiff = Jdiff_interp*gridReachAvgSVF'; %diffuse in canyon
        
        Jsn_total = Jreflected+Jdir+Jdiff; % total Jsn at water surface
        
        radiationFactor = 1-(Jsn_total./Jsn_interp);
        radiationFactor(isnan(radiationFactor))=1;
        era5.selected.(eraFields{i}).RadiationFactors = radiationFactor;
        
        %% Write out to file - radiation factor
        cd('output_data')
        filename = char(strcat('era5_radiationFactors_',eraFields(i),'.csv'));
        fid = fopen(filename,'W');
        fprintf(fid,'%s,','Datetime');
        for j = 1:length(gridRiverKMs)
            fprintf(fid,'L_%s,',num2str(gridRiverKMs(j)));
        end
        
        fprintf(fid,'\n');
        for j = 1:length(writeDateTimesStr) % loop for datetimes
            % print datetime
            fprintf(fid,'%s,',writeDateTimesStr(j,:));
            for jj = 1:size(radiationFactor,2) % loop for averaged incidence
                % print shade factor for each element
                fprintf(fid,'%1.2g,',radiationFactor(j,jj));
            end
            fprintf(fid,'\n');
            
            if mod(j,round(size(radiationFactor,1)/100))==0 %write out percent progress to command window
                percentComplete = round(j/size(radiationFactor,1)*100);
                fprintf('~%d%% complete\n',percentComplete);
            end
        end
        fclose(fid);
        cd(mfileDir);
    end
    
    cd(mfileDir)
    save(matlabDataFile,'era5','uInd','subs')
else
    load('matlabDataFile.mat')
end

%% Check net shortwave radiation
% for i = 24% 1:length(eraFields)
%     netJsn = era5.selected.(eraFields{i}).RadiationFactors .* era5.selected.(eraFields{i}).SolarRadiation{2};
%     dt = era5.selected.(eraFields{i}).SolarRadiation{1};
%     dateRange = [datetime('11/19/2008'):hours(1):datetime('11/25/2008')];
%     ind = ismember(dt,dateRange);
%     netJsn = netJsn(ind,:);
%     dt = dt(ind);
%     x = 1:size(netJsn,2);
%     
%     figure()
%     hold on
%     for j = 1:length(x)
%         plot(dt,netJsn(:,j))
%     end
%     ax = gca;
%     ax.Title.String = eraFields{i};
%     close
%     close
%     
%     %     contourf(x,datenum(dt),netJsn)
%     %     ax = gca;
%     %     c = ax.Children;
%     %     c.LineStyle = 'none';
%     %     ax.XAxis.TickDirection = 'out';
%     %     ax.YLabel.String = 'Date';
%     %     ax.YTickLabel = datestr(ax.YTick,'mm/DD/YYYY');
%     %     ax.YAxis.TickDirection = 'out';
%     %     ax.Title.String = eraFields{i};
%     %     cb = colorbar;
%     %     cb.Label.String = 'Net Jsn';
%     %     cb.FontSize = ax.YLabel.FontSize;
%     %     pause();
%     %     close
% end

