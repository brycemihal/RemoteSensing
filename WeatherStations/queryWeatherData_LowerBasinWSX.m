%% Connect to the database

dbPath = 'C:\Users\Echo\Neilson Research Dropbox\Projects\ColoradoRiverBasinProject\Data_mgmt\Databases\';
dbName = 'LowerCORiverBasin_full.db'; % change if working in the Upper Basin
conn = setSQLitedriver(dbName,dbPath);

%% query site data

qstr = "SELECT * FROM sites WHERE Comments = 'LowerBasin' AND SiteType = 'Climate'";
curs = exec(conn,char(qstr));
curs = fetch(curs);
siteData = curs.Data;

siteID = cell2mat(siteData(:,1));

%% get air temperature data (vID = 1)
for i = 1:length(siteID)
qstr = strcat({'SELECT LocalDateTime,DataValue FROM datavalues '},...
    {'WHERE VariableID = 1 and SiteID = $siteid$ and DataValue <> -9999'});

qstr = strrep(qstr,'$siteid$',num2str(siteID((i))));

curs = exec(conn,char(qstr));
curs = fetch(curs);
atData{i} = curs.Data;
end

%% get relative humidity data (vID = 3)
for i = 1:length(siteID)
qstr = strcat({'SELECT LocalDateTime,DataValue FROM datavalues '},...
    {'WHERE VariableID = 3 and SiteID = $siteid$ and DataValue <> -9999'});

qstr = strrep(qstr,'$siteid$',num2str(siteID((i))));

curs = exec(conn,char(qstr));
curs = fetch(curs);
rhData{i} = curs.Data;
end

%% get solar radiation data (vID = 4)
for i = 1:length(siteID)
qstr = strcat({'SELECT LocalDateTime,DataValue FROM datavalues '},...
    {'WHERE VariableID = 4 and SiteID = $siteid$ and DataValue <> -9999'});

qstr = strrep(qstr,'$siteid$',num2str(siteID((i))));

curs = exec(conn,char(qstr));
curs = fetch(curs);
srData{i} = curs.Data;
end

%% get wind speed data (vID = 5)
for i = 1:length(siteID)
qstr = strcat({'SELECT LocalDateTime,DataValue FROM datavalues '},...
    {'WHERE VariableID = 5 and SiteID = $siteid$ and DataValue <> -9999'});

qstr = strrep(qstr,'$siteid$',num2str(siteID((i))));

curs = exec(conn,char(qstr));
curs = fetch(curs);
wsData{i} = curs.Data;
end

%% Build structure 
siteNames = strrep(strrep(siteData(:,3),'-',''),{' '},'_');

siteNames{20} = 'Mile_n10_Ferry_Swale';
siteNames{21} = 'Mile_0_5_Lees_Ferry';
siteNames{22} = 'Mile_24_5_Upper';
siteNames{23} = 'Mile_24_5_Lower';
siteNames{24} = 'Mile_58_Upper';
siteNames{25} = 'Mile_58_Lower';
siteNames{26} = 'Mile_60';
siteNames{28} = 'Mile_70_Upper';
siteNames{29} = 'Mile_70_Lower';
siteNames{30} = 'Mile_125_5_Fossil';
siteNames{32} = 'Mile_203';
siteNames{33} = 'Mile_223';

siteAttributes = {'SiteID','SiteCode','Latitude','Longitude','State',...
    'Elevation','Data'};

for i = 1:length(siteNames)
    lowerBasinWSX.(siteNames{i}).(siteAttributes{1}) = cell2mat(siteData(i,1));
    lowerBasinWSX.(siteNames{i}).(siteAttributes{2}) = char(siteData(i,2));
    lowerBasinWSX.(siteNames{i}).(siteAttributes{3}) = cell2mat(siteData(i,4));
    lowerBasinWSX.(siteNames{i}).(siteAttributes{4}) = cell2mat(siteData(i,5));
    lowerBasinWSX.(siteNames{i}).(siteAttributes{5}) = char(siteData(i,13));
    lowerBasinWSX.(siteNames{i}).(siteAttributes{6}) = cell2mat(siteData(i,7));
    lowerBasinWSX.(siteNames{i}).(siteAttributes{7}).AirTemp = atData{i};
    lowerBasinWSX.(siteNames{i}).(siteAttributes{7}).RelativeHumidity = rhData{i};
    lowerBasinWSX.(siteNames{i}).(siteAttributes{7}).WindSpeed = wsData{i};
    lowerBasinWSX.(siteNames{i}).(siteAttributes{7}).SolarRadiation = srData{i};
end
