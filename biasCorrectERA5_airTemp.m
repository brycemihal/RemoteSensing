%% Script to calculate IDW factors for ERA5 Grids
% Bryce Mihalevich
clear; clc; close all;

% Relative path to files
n = mfilename;
p = mfilename('fullpath');
cd(p(1:end-length(n)))
mfileDir = pwd;

%% Import data
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\wsx_comparison\era5Land_wsx_noPage.mat');
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\grids_to_model_elements\era5Land.mat')
load('C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Weather Data\Reanalysis\ERA5Land\bias_correction\era5_airTempBiasCorrectionFactors_noPage.mat')

%% Set up data
% hand pick weather stations to use (GCMRC sites)
gcmrc_sites = [2 3 8 12 13 17 18 19 20 24];
uInd_wsx = uInd_wsx(gcmrc_sites);

%get all unique pairs of river grid points
[era5_lat,era5_lon] = meshgrid(era5.latitude,era5.longitude);
era5_lat = era5_lat(:);
era5_lon = era5_lon(:);

%get all unique pairs of weather grid points
[era5_wsx_lat,era5_wsx_lon] = meshgrid(era5_wsx.latitude,era5_wsx.longitude);
era5_wsx_lat = era5_wsx_lat(:);
era5_wsx_lon = era5_wsx_lon(:);

%% Perform Bias Correction
% loop for river grids
for i = 1:length(uInd)
    % calculate distance between era river grid and weather station grid
    for j = 1:length(uInd_wsx)
        % calculate ellipsoid distance
        d(j,1) = geoddistance(era5_lat(uInd(i)),era5_lon(uInd(i)),...
            era5_wsx_lat(uInd_wsx(j)),era5_wsx_lon(uInd_wsx(j)));
    end
    
    % if river grid includes a weather site don't weight
    if min(d) == 0
        weightFactor = double(d==0);
    else
        % calculate weight for each wsx
        n = 3; %weighting factor. Increase to 2 to give higher wieght to nearer stations
        wi = 1./(d.^n);
        weightFactor = (wi./sum(wi)); % each row contains weight factors for each gcmrc site
    end
    
    % apply bias correction to grid cell
    era5 = applyBiasCorrection(i,era5,weightFactor,linearOffset,varianceMultiplyer);
end

%% Save output

save('era5_biasCorrected.mat','era5','uInd','subs')

%% Write selected ERA5 grids to files
gridID = fieldnames(era5.biasCorrected);
params = fieldnames(era5.biasCorrected.(gridID{1}));
pdatetime = datestr(era5.biasCorrected.(gridID{1}).(params{1}){1},'mm/dd/yyyy HH:MM:SS');

save2Dir = fullfile(mfileDir,'bias_corrected_output_for_model');
if ~exist(save2Dir,'dir')
    mkdir(save2Dir);
end
cd(save2Dir);

% % model element to grid reference file
% metNames = {'AIR_TEMPERATURE','WIND_SPEED','','RELATIVE_HUMIDITY'}; %no name for Solar Radiation in input file
% refFile = strrep(riverFile,'.csv','_ERA5.txt');
% fid = fopen(refFile,'w');
% fprintf(fid,'[METEOROLOGY]\n');
% fprintf(fid,';;START_ELEMENT    END_ELEMENT    VARIABLE           TYPE         TIMESERIES\n');
% fprintf(fid,';;==========================================================================\n');
% fclose(fid);
% fid = fopen(refFile,'a');
% for p = 1:length(params)
%     for i = 1:length(subs)
%         gID = gridID{subs(i)};
%         seriesPath = strcat('./ERA5/ERA5_',params{p},'_',gID,'_bc.csv');
%         seriesName = strrep(strrep(seriesPath,'.csv',''),'./ERA5/','');
%         printstr = strcat('L_',num2str(rivCentroidKM(i)),{'  '},...
%             'L_',num2str(rivCentroidKM(i)),{'  '},...
%             metNames{p},{'  '},'TIMESERIES',{'  '},seriesName);
%         fprintf(fid,'%s\n',char(printstr));
%     end
%     fprintf(fid,';;==========================================================================\n');
% end
% fclose(fid);
% 
% % model element to grid reference file
% refFile = strrep(riverFile,'.csv','_ERA5_timeseries.txt');
% fid = fopen(refFile,'w');
% fprintf(fid,'[TIMESERIES]\n');
% fprintf(fid,';;NAME                              FILEPATH\n');
% fprintf(fid,';;==========================================================================\n');
% fclose(fid);
% fid = fopen(refFile,'a');
% for p = 1:length(params)
%     for i = 1:length(gridID)
%         seriesPath = strcat('./ERA5/ERA5_',params{p},'_',gridID{i},'_bc.csv');
%         seriesName = strrep(strrep(seriesPath,'.csv',''),'./ERA5/','');
%         printstr = strcat(seriesName,{'  '},seriesPath);
%         fprintf(fid,'%s\n',char(printstr));
%     end
%     fprintf(fid,';;==========================================================================\n');
% end
% fclose(fid);

%% write parameters for each unique grid
for i = 1:length(gridID)
    for p = 1:length(params)
        filename = strcat('ERA5_',params{p},'_',gridID{i},'_bc.csv');
        fid = fopen(filename,'w');
        fprintf(fid,strcat('Datetime (MST),',filename)); %limit the length of the filename, seems to crash swmm when too long
        fprintf(fid,'\n');
        fclose(fid);
        
        pdata = round(era5.biasCorrected.(gridID{i}).(params{p}){2},2);
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

%% Functions:
function[era5] = applyBiasCorrection(i,era5,weightFactor,linearOffset,varianceMultiplyer)
gridID = fieldnames(era5.selected);
dt = era5.selected.(gridID{i}).AirTemp{1}; %datetimes
dv = era5.selected.(gridID{i}).AirTemp{2}; %datavalues

% weighted linear correction
linCorrFactor = sum(linearOffset.*weightFactor',2);

% weighted variance multiplyer
varCorrFactor = sum(varianceMultiplyer.*weightFactor',2);

% apply linear offset correction
linearCorrection = dv + linCorrFactor(dt.Month);

% apply variance scaling correction
linearCorrectionMonthly = accumarray(dt.Month,linearCorrection,[],@mean);
varianceCorrection = ((linearCorrection - linearCorrectionMonthly(dt.Month))...
    .* varCorrFactor(dt.Month))+ linearCorrectionMonthly(dt.Month);

era5.biasCorrected.(gridID{i}).AirTemp = {dt varianceCorrection};
end


