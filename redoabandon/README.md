
# Time series of Ocean Color data

Welcome to my script


I would like to introduce you to how to create a time series of ocean color data.

>>Before anything, go to NASA's ocean color database (HERE: https://oceandata.sci.gsfc.nasa.gov/api/file_search/) and download the netCDF files for each parameter you are interested in and for every 8day, monthly, or annual time interval (your choice). HDF files for global matrix (OC_coverageExample.mat), net primary productivity, and mixed layer depth are also available (HERE: http://orca.science.oregonstate.edu/)

>> You will want to organize files into folders by parameters within your main directory.

>> Coding process:
First read in latitude and longitude coordinates:

```matlab
%% Load Data and Define Time Range
% Load lon and lat data from the .mat file
x = load('OC_coverageExample.mat');

```

>>Then will define each collection of parameters by loading in the folders organized in the second step

```matlab
% Load all files in the directories
filechl    = dir('chl/AQUA*chlor_a.9km.nc');
filepar    = dir('par/AQUA*.par.9km.nc');
filebbp443 = dir('bbp443/AQUA*.9km.nc');
filesst    = dir('sst/AQUA*.9km.nc');
filemld    = dir('MLD/mld*.hdf');
filenpp    = dir('NPP/npp*.hdf');

```

>>Define your time series as shown below
```matlab
%% Define Time Range

startDate  = datetime(2002,07,04); %customize to your analysis
endDate    = datetime(2021,12,31); %customize to your analysis
timeStep   = caldays(8); %customize to your analysis
timeVector = startDate:timeStep:endDate;
num_times  = length(timeVector);

% Convert the datetime vector to numeric values for polyfit
time_num = datenum(timeVector)';
```

>>Ensure your coordinates align with the OC_coverageExample.mat file and define your region of interest
```matlab
%% Load lon/lat and Set Coordinate Ranges
lon = ncread("chl/AQUA_MODIS.20020704_20020711.L3m.8D.CHL.chlor_a.9km.nc", 'lon');
lat = ncread("chl/AQUA_MODIS.20020704_20020711.L3m.8D.CHL.chlor_a.9km.nc", 'lat');

% Define coordinate ranges
lon_ranges = [-150, -140]; %customize to your analysis
lat_ranges = [50, 55]; %customize to your analysis

% Get indices for coordinate ranges
lon_idx = find(lon >= lon_ranges(1) & lon <= lon_ranges(2));
lat_idx = find(lat >= lat_ranges(1) & lat <= lat_ranges(2));
```

>>Then create empty arrays for each parameter, you will need fill these with the mean/median/min/max values calculated later
```matlab
%% Preallocate arrays for mean values
avg_chla_values   = NaN(num_times, 1);
avg_par_values    = NaN(num_times, 1);
avg_sst_values    = NaN(num_times, 1);
avg_bbp443_values = NaN(num_times, 1);
avg_cphyto_values = NaN(num_times, 1);
avg_mld_values    = NaN(num_times, 1);
avg_npp_values    = NaN(num_times, 1);
```

>>Now you are set up to write for loops or functions to calculate values for each parameter, here I used for loops and the tic toc function in matlab to time each section of the script
```matlab
%% Process Chlorophyll Data

tic
for t = 1:num_times
    if t <= length(filechl)
        chlor_a = ncread(fullfile('chl', filechl(t).name), 'chlor_a');
        region_chla = chlor_a(lon_idx, lat_idx);
        avg_chla_values(t) = mean(region_chla(:), 'omitnan');
    end
end
toc

%% Process NPP Data

tic
for t = 1:num_times
    if t <= length(filenpp)
        npp = hdfread(fullfile('NPP', filenpp(t).name), 'npp');
        region_npp = npp(lon_idx, lat_idx);
        avg_npp_values(t) = mean(region_npp(:), 'omitnan');
    end
end
toc

%% Process PAR Data

tic
for t = 1:num_times
    if t <= length(filepar)
        par = ncread(fullfile('par', filepar(t).name), 'par');
        region_par = par(lon_idx, lat_idx);
        avg_par_values(t) = mean(region_par(:), 'omitnan');
    end
end
toc
%% Process bbp443 and Phytoplankton Biomass (Cphyto) Data

tic
for t = 1:num_times
    if t <= length(filebbp443)
        bbp443 = ncread(fullfile('bbp443', filebbp443(t).name), 'bbp_443');
        Cphyto = 13000 .* bbp443 + 4.55; % Calculate Cphyto
        region_bbp443 = bbp443(lon_idx, lat_idx);
        region_cphyto = Cphyto(lon_idx, lat_idx);
        avg_bbp443_values(t) = mean(region_bbp443(:), 'omitnan');
        avg_cphyto_values(t) = mean(region_cphyto(:), 'omitnan');
    end
end
toc

%% Process Sea Surface Temperature (SST) Data

tic
for t = 1:num_times
    if t <= length(filesst)
        sst = ncread(fullfile('sst', filesst(t).name), 'sst');
        region_sst = sst(lon_idx, lat_idx);
        avg_sst_values(t) = mean(region_sst(:), 'omitnan');
    end
end
toc

%% Process Mixed Layer Depth (MLD) Data

tic
for t = 1:num_times
    if t <= length(filemld)
        mld = hdfread(fullfile('MLD', filemld(t).name), 'mld');
        mld = mld';
        region_mld = mld(lon_idx, lat_idx);
        avg_mld_values(t) = mean(region_mld(:), 'omitnan');
    end
end

% Convert the mean MLD to negative values
avg_mld_values = -avg_mld_values;
toc

```
>> Now that you have time series variables in your workspace, you can move on to creating figures and finding some statistics. This is what I did but there are many different ways to customize your figures, explore.
```matlab
%% Plotting and Saving Chlorophyll Data with Best Fit Line

tic
figure;
plot(timeVector, avg_chla_values, 'Color', [0.1 0.6470 0.7410], 'DisplayName', 'Chlorophyll Data');
hold on;
% Compute best fit line
p_chla = polyfit(time_num, avg_chla_values, 1);
chla_fit = polyval(p_chla, time_num);
plot(timeVector, chla_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');
title('Gulf Alask Average Chlorophyll Concentration Over Time');
xlabel('Date');
ylabel('Average Chlorophyll (mg chl^-^m^3)');

% Add callouts for max and min points
[~, maxIdx] = max(avg_chla_values);
[~, minIdx] = min(avg_chla_values);
plot(timeVector(maxIdx), avg_chla_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_chla_values(maxIdx), sprintf('Max: %.2f', avg_chla_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_chla_values(minIdx), 'ro');
text(timeVector(minIdx), avg_chla_values(minIdx), sprintf('Min: %.2f', avg_chla_values(minIdx)), 'VerticalAlignment', 'top');
legend show;
hold off;
saveas(gcf, 'np_average_chlorophyll.png');
toc

%% Plotting and Saving Phytoplankton Biomass (Cphyto) Data with Best Fit Line

tic
figure;
plot(timeVector, avg_cphyto_values, 'Color', [0.6 0.1 0.8], 'DisplayName', 'PhytoBiomass Data');
hold on;
% Compute best fit line
p_cphyto = polyfit(time_num, avg_cphyto_values, 1);
cphyto_fit = polyval(p_cphyto, time_num);
plot(timeVector, cphyto_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');
title('Gulf Alaska Average PhytoBiomass Concentration Over Time');
xlabel('Date');
ylabel('Average PhytoBiomass (mg C^-^L)');

% Add callouts for max and min points
[~, maxIdx] = max(avg_cphyto_values);
[~, minIdx] = min(avg_cphyto_values);
plot(timeVector(maxIdx), avg_cphyto_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_cphyto_values(maxIdx), sprintf('Max: %.2f', avg_cphyto_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_cphyto_values(minIdx), 'ro');
text(timeVector(minIdx), avg_cphyto_values(minIdx), sprintf('Min: %.2f', avg_cphyto_values(minIdx)), 'VerticalAlignment', 'top');
legend show;
hold off;
saveas(gcf, 'np_average_phytobiomass.png');
toc

%% Plotting and Saving PAR Data with Best Fit Line

tic
figure;
plot(timeVector, avg_par_values, 'Color', [0.8500 0.4250 0.0980], 'DisplayName', 'PAR Data');
hold on;
% Compute best fit line
p_par = polyfit(time_num, avg_par_values, 1);
par_fit = polyval(p_par, time_num);
plot(timeVector, par_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');
title('Gulf Alaska Average PAR Over Time');
xlabel('Date');
ylabel('Average PAR (mol photons^-^m^2^-^d)');
toc

% Add callouts for max and min points
[~, maxIdx] = max(avg_par_values);
[~, minIdx] = min(avg_par_values);
plot(timeVector(maxIdx), avg_par_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_par_values(maxIdx), sprintf('Max: %.2f', avg_par_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_par_values(minIdx), 'ro');
text(timeVector(minIdx), avg_par_values(minIdx), sprintf('Min: %.2f', avg_par_values(minIdx)), 'VerticalAlignment', 'top');
legend show;
hold off;
saveas(gcf, 'np_average_par.png');
toc

%% Plotting and Saving SST Data with Best Fit Line

tic
figure;
plot(timeVector, avg_sst_values, 'Color', [0.8350 0.0780 0.2840], 'DisplayName', 'SST Data');
hold on;
% Compute best fit line
p_sst = polyfit(time_num, avg_sst_values, 1);
sst_fit = polyval(p_sst, time_num);
plot(timeVector, sst_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');
title('Gulf Alaska Average Sea Surface Temperature Over Time');
xlabel('Date');
ylabel('Average SST (°C)');

% Add callouts for max and min points
[~, maxIdx] = max(avg_sst_values);
[~, minIdx] = min(avg_sst_values);
plot(timeVector(maxIdx), avg_sst_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_sst_values(maxIdx), sprintf('Max: %.2f°C', avg_sst_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_sst_values(minIdx), 'ro');
text(timeVector(minIdx), avg_sst_values(minIdx), sprintf('Min: %.2f°C', avg_sst_values(minIdx)), 'VerticalAlignment', 'top');
legend show;
hold off;
saveas(gcf, 'np_average_sst.png');
toc

%% Plotting and Saving MLD Data with Best Fit Line

tic
figure;
plot(timeVector, avg_mld_values', 'Color', [0.3 0.8 0.3], 'DisplayName', 'MLD Data');
hold on;

% Compute best fit line
p_mld = polyfit(time_num, avg_mld_values, 20);
mld_fit = polyval(p_mld, time_num);
plot(timeVector, mld_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');
title('Gulf Alaska Average Mixed Layer Depth Over Time');
xlabel('Date');
ylabel('Average MLD (m)');

% For MLD the "max" value is the most negative, so we use min() to call it out accordingly.
[~, maxIdx] = min(avg_mld_values);
[~, minIdx] = max(avg_mld_values);
plot(timeVector(maxIdx), avg_mld_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_mld_values(maxIdx), sprintf('Max: %.2f m', avg_mld_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_mld_values(minIdx), 'ro');
text(timeVector(minIdx), avg_mld_values(minIdx), sprintf('Min: %.2f m', avg_mld_values(minIdx)), 'VerticalAlignment', 'top');
legend show;
hold off;
saveas(gcf, 'np_average_mld.png');
toc

%% Plotting and Saving NPP Data with Best Fit Line

tic
figure;
plot(timeVector, avg_npp_values, 'Color', [0.7 0.1 0.4], 'DisplayName', 'NPP Data');
hold on;

% Compute best fit line
p_npp = '(time_num; avg_npp_values; 1)';
npp_fit = polyval(p_npp, time_num);
plot(timeVector, npp_fit, 'k--', 'LineWidth', 2, 'DisplayName', 'Best Fit Line');
title('Tit Average Net Primary Production Over Time');
xlabel('Date');
ylabel('Average NPP  (g Cm^-^2d^-^1)');

% For NPP the "max" value is the most negative, so we use min() to call it out accordingly.
[~, maxIdx] = min(avg_npp_values);
[~, minIdx] = max(avg_npp_values);
plot(timeVector(maxIdx), avg_npp_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_npp_values(maxIdx), sprintf('Max: %.2f m', avg_npp_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_npp_values(minIdx), 'ro');
text(timeVector(minIdx), avg_npp_values(minIdx), sprintf('Min: %.2f m', avg_npp_values(minIdx)), 'VerticalAlignment', 'top');
legend show;
hold off;
saveas(gcf, 'np_average_npp.png');
toc

```
>>If you desire to work in excel, or other platform, create a table of all the dataset created in this code
```matlab
%% Create a Table and Save as Excel File

tic
T = table(timeVector', avg_chla_values, avg_par_values, avg_sst_values, avg_bbp443_values, avg_cphyto_values, avg_mld_values, avg_npp_values ...
    ,'VariableNames', {'Time', 'Average_Chlorophyll', 'Average_PAR', 'Average_SST','Average_BBP443', 'Average_Cphyto', 'Average_MLD', 'Average_NPP'});
    writetable(T, 'np_average_values_over_time.xlsx'); %you can make this a csv file for platform versatility.
toc

```
