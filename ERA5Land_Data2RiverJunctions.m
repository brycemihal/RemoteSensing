%% Script to generate csv files for each grid that intersects the river
% Bryce Mihalevich
% 1/24/2020

% info
% https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview

% Relative path to files
n = mfilename;
p = mfilename('fullpath');
cd(p(1:end-length(n)))
mfileDir = pwd;
% cd ..
% era5Dir = pwd;
%% Import river element junction data
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

%% Import data from the netcdf
matlabDataFile = 'era5Land2.mat';
if ~exist(matlabDataFile,'file')>0  
    cd ..
    cd('raw_data')
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
        
        %% Determine what reanalysis grids correspond to river elements
        %get all unique pairs of grid points
        [era5_lat,era5_lon] = meshgrid(era5.latitude,era5.longitude);
        era5_lat = era5_lat(:);
        era5_lon = era5_lon(:);
        
        % convert UTM to lat lon
        [rivCentroid_lat, rivCentroid_lon] = utmups_inv(rivCentroidX,rivCentroidY, 12, 1);
                
        for j = 1:length(rivCentroidKM)
            % caclulate distance between river element and RS cells and grab index of minimum distances
            [~,ind(j)] = min(sqrt((rivCentroid_lon(j)-era5_lon).^2 + (rivCentroid_lat(j)-era5_lat).^2));
        end
        
        [uInd,~,subs] = unique(ind);
        
        %% Plot visualization
%         close all
%         figure()
%         hold on        
%         % add era5 grid boxes
%         lat_spacing = mode(diff(unique(era5_lat)))/2;
%         lon_spacing = mode(diff(unique(era5_lon)))/2;
%         lat_lines = [min(era5_lat)-lat_spacing; unique(era5_lat)+lat_spacing];
%         lon_lines = [min(era5_lon)-lon_spacing; unique(era5_lon)+lon_spacing];
%         for j = 1:length(lat_lines)
%             plot([min(lon_lines) max(lon_lines)],[lat_lines(j) lat_lines(j)],'k')
%         end
%         for j = 1:length(lon_lines)
%             plot([lon_lines(j),lon_lines(j)],[min(lat_lines) max(lat_lines)],'k')
%         end
%         % add patch for selected boxes
%         for j = 1:length(uInd)
%             xq = [era5_lon(uInd(j))-lon_spacing era5_lon(uInd(j))+lon_spacing era5_lon(uInd(j))+lon_spacing era5_lon(uInd(j))-lon_spacing]; 
%             yq = [era5_lat(uInd(j))-lat_spacing era5_lat(uInd(j))-lat_spacing era5_lat(uInd(j))+lat_spacing era5_lat(uInd(j))+lat_spacing];
%             p = patch(xq,yq,'red');
%             p.FaceColor = [1 0 0];
%             p.FaceAlpha = .5;
%             t = text();
%             t.String = num2str(j);
%             t.Units = 'data';
%             t.FontSize = 10;
%             t.FontWeight = 'bold';
%             t.Position = [era5_lon(uInd(j)) era5_lat(uInd(j)) 0];
%             t.HorizontalAlignment = 'center';
%             t.VerticalAlignment = 'middle';
%             t.Clipping = 'on';
%         end
% %         plot(era5_lon,era5_lat,'ok')
% %         plot(rivCentroid_lon,rivCentroid_lat,'oy','MarkerSize',2)       
% %         plot(site_lon,site_lat,'ok','MarkerFaceColor','k')
%         plot_google_map('MapType','terrain','ShowLabels',0)
%         ax = gca;
%         ax.YLabel.String = 'Latitude';
%         ax.XLabel.String = 'Longitude';
%         
% %         str = ['grid = ERA dataset' newline 'black dots = weather stations' newline 'red boxes = nearest ERA grid'];
% %         t = text();
% %         t.String = str;
% %         t.Units = 'normalized';
% %         t.Position = [0.01 .9 0];
% %         t.FontSize = 12;
% %         t.BackgroundColor = 'w';
        
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
            % note: sometimes there are NaN values.
            qSSRD = squeeze(cell2mat({era5.solar_down(era5_lon_ind,era5_lat_ind,:)}))/3600;
            
            % loop through each timestep
            qSolar(1,1) = 0;
            for jj = 2:length(era5.time_utc) %qSSRD is the accumulation of Jsn through a day. 
                if hour(era5.time_utc(jj))==1 %need to get Jsn per hour
                    qSolar(jj,1) = qSSRD(jj);
                else
                    qSolar(jj,1) = qSSRD(jj)- qSSRD(jj-1);
                end
            end
            
            if sum(isnan(qSolar))>0
                qSolar(isnan(qSolar))=0; % set NaN values to zero. Have only seen this at night time. 
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
    
    %% Check for NaN values
    gridID = fieldnames(era5.selected);
    params = fieldnames(era5.selected.(gridID{1}));

    nanStop = false;
    for i = 1:length(gridID)
        for j = 1:length(params)
            dt = era5.selected.(gridID{i}).(params{j}){1};
            nanVals = isnan(era5.selected.(gridID{i}).(params{j}){2});
            if sum(nanVals) > 0            
                warning(sprintf('NaN value found in %s of %s',params{j},gridID{i})) 
                disp(dt(nanVals))
                nanStop = true;
            end
        end
    end
    if nanStop
        error('Fix NaN values')
    end
    
    cd(mfileDir)
    save(matlabDataFile,'era5','uInd','subs')
else
    load(matlabDataFile)
end

%% Write selected ERA5 grids to files
gridID = fieldnames(era5.selected);
params = fieldnames(era5.selected.(gridID{1}));
pdatetime = datestr(era5.selected.(gridID{1}).(params{1}){1},'mm/dd/yyyy HH:MM:SS');

save2Dir = fullfile(mfileDir,'output_for_model');
if ~exist(save2Dir,'dir')
    mkdir(save2Dir);
end
cd(save2Dir);

% model element to grid reference file
metNames = {'AIR_TEMPERATURE','WIND_SPEED','','RELATIVE_HUMIDITY'}; %no name for Solar Radiation in input file
refFile = strrep(riverFile,'.csv','_ERA5.txt');
fid = fopen(refFile,'w');
fprintf(fid,'[METEOROLOGY]\n');
fprintf(fid,';;START_ELEMENT    END_ELEMENT    VARIABLE           TYPE         TIMESERIES\n');
fprintf(fid,';;==========================================================================\n');
fclose(fid);
fid = fopen(refFile,'a');
for p = 1:length(params)
    for i = 1:length(subs)
        gID = gridID{subs(i)};
        seriesPath = strcat('./ERA5/ERA5_',params{p},'_',gID,'.csv');
        seriesName = strrep(strrep(seriesPath,'.csv',''),'./ERA5/','');
        printstr = strcat('L_',num2str(rivCentroidKM(i)),{'  '},...
            'L_',num2str(rivCentroidKM(i)),{'  '},...
            metNames{p},{'  '},'TIMESERIES',{'  '},seriesName);
        fprintf(fid,'%s\n',char(printstr));
    end
    fprintf(fid,';;==========================================================================\n');
end
fclose(fid);

% model element to grid reference file
refFile = strrep(riverFile,'.csv','_ERA5_timeseries.txt');
fid = fopen(refFile,'w');
fprintf(fid,'[TIMESERIES]\n');
fprintf(fid,';;NAME                              FILEPATH\n');
fprintf(fid,';;==========================================================================\n');
fclose(fid);
fid = fopen(refFile,'a');
for p = 1:length(params)
    for i = 1:length(gridID)
        seriesPath = strcat('./ERA5/ERA5_',params{p},'_',gridID{i},'.csv');
        seriesName = strrep(strrep(seriesPath,'.csv',''),'./ERA5/','');
        printstr = strcat(seriesName,{'  '},seriesPath);
        fprintf(fid,'%s\n',char(printstr));
    end
    fprintf(fid,';;==========================================================================\n');
end
fclose(fid);

% write parameters for each unique grid
for i = 1:length(gridID)
    for p = 3%1:length(params)
        filename = strcat('ERA5_',params{p},'_',gridID{i},'.csv');
        fid = fopen(filename,'w');
        fprintf(fid,strcat('Datetime (MST),',filename)); %limit the length of the filename, seems to crash swmm when too long
        fprintf(fid,'\n');
        fclose(fid);
        
        pdata = round(era5.selected.(gridID{i}).(params{p}){2},2);
        fid = fopen(filename,'A');
        for j = 1:length(pdata)
            fprintf(fid,'%s,%g\n',pdatetime(j,:),pdata(j));
            
            if mod(j,round(length(pdata)/25))==0 %write out percent progress to command window
                percentComplete = round(j/length(pdata)*100);
                pp = sprintf('~%d%% complete',percentComplete);
                disp(pp)
            end
        end
        fclose(fid);
    end
end

cd ..
