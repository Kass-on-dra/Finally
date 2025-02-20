 clc; clearvars;
%%
% Load lon and lat data from the .mat file
x = load('OC_coverageExample.mat');

% Load all files in the directories
filechl = dir('chl/AQUA*chlor_a.9km.nc');
filepar = dir('par/AQUA*.par.9km.nc');
filebbp443 = dir('bbp443/AQUA*.9km.nc');
filesst = dir('sst/AQUA*.9km.nc');
filemld = dir('MLD/mld*.hdf');

% Use time range and calculate expected number of files
startDate = datetime(2002, 07, 04);
endDate = datetime(2021, 12, 31);
timeStep = caldays(8);
timeVector = startDate:timeStep:endDate;
num_times = length(timeVector);

% Load lon/lat & mesh
lon = ncread("chl/AQUA_MODIS.20020704_20020711.L3m.8D.CHL.chlor_a.9km.nc", 'lon');
lat = ncread("chl/AQUA_MODIS.20020704_20020711.L3m.8D.CHL.chlor_a.9km.nc", 'lat');

% Define coordinate ranges
lon_ranges = [-60, -50];
lat_ranges = [40, 45];

% Get indices for coordinate ranges
lon_idx = find(lon >= lon_ranges(1) & lon <= lon_ranges(2));
lat_idx = find(lat >= lat_ranges(1) & lat <= lat_ranges(2));
%% Preallocate for average values
avg_chla_values = NaN(num_times, 1);
avg_par_values = NaN(num_times, 1);
avg_sst_values = NaN(num_times, 1);
avg_bbp443_values = NaN(num_times, 1);
avg_cphyto_values = NaN(num_times, 1);
avg_mld_values = NaN(num_times, 1);

%% Process Chlorophyll data
for t = 1:num_times
    if t <= length(filechl)
        chlor_a = ncread(fullfile('chl', filechl(t).name), 'chlor_a');
        region_chla = chlor_a(lon_idx, lat_idx);
        avg_chla_values(t) = mean(region_chla(:), 'omitnan');

    end
end

%% Process PAR data
for t = 1:num_times
    if t <= length(filepar)
        par = ncread(fullfile('par', filepar(t).name), 'par');
        region_par = par(lon_idx, lat_idx);
        avg_par_values(t) = mean(region_par(:), 'omitnan');
    end
end

%% Process bbp443 and Cphyto data
for t = 1:num_times
    if t <= length(filebbp443)
        bbp443 = ncread(fullfile('bbp443', filebbp443(t).name), 'bbp_443');
        Cphyto = 12128 .* bbp443 + 0.59; % Calculate Cphyto
        region_bbp443 = bbp443(lon_idx, lat_idx);
        region_cphyto = Cphyto(lon_idx, lat_idx);
        avg_bbp443_values(t) = mean(region_bbp443(:), 'omitnan');
        avg_cphyto_values(t) = mean(region_cphyto(:), 'omitnan');
    end
end

%% Process SST data
for t = 1:num_times
    if t <= length(filesst)
        sst = ncread(fullfile('sst', filesst(t).name), 'sst');
        region_sst = sst(lon_idx, lat_idx);
        avg_sst_values(t) = mean(region_sst(:), 'omitnan');
    end
end

%% Process MLD data
for t = 1:num_times
    if t <= length(filemld)
        mld = hdfread(fullfile('MLD', filemld(t).name), 'mld');
        mld = mld';
        region_mld = mld(lon_idx, lat_idx);
        avg_mld_values(t) = mean(region_mld(:), 'omitmissing');
    end
end

% Convert the mean MLD to negative values
avg_mld_values = -avg_mld_values;

%% Plotting and Saving Chlorophyll
figure;
plot(timeVector, avg_chla_values, 'Color', [0.1 0.6470 0.7410]);
title('Tit Average Chlorophyll Concentration Over Time');
xlabel('Date');
ylabel('Average Chlorophyll (mg chl^-^m^3)');

% Add callouts
[~, maxIdx] = max(avg_chla_values);
[~, minIdx] = min(avg_chla_values);
hold on;
plot(timeVector(maxIdx), avg_chla_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_chla_values(maxIdx), ...
    sprintf('Max: %.2f', avg_chla_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_chla_values(minIdx), 'ro');
text(timeVector(minIdx), avg_chla_values(minIdx), ...
    sprintf('Min: %.2f', avg_chla_values(minIdx)), 'VerticalAlignment', 'top');
hold off;

% Save figure
saveas(gcf, 'tit_average_chlorophyll.png');

%% Plotting and Saving Phytoplankton Biomass
figure;
plot(timeVector, avg_cphyto_values, 'Color', [0.6 0.1 0.8]);
title('Tit Average PhytoBiomass Concentration Over Time');
xlabel('Date');
ylabel('Average PhytoBiomass (mg C^-^L^m)');

% Add callouts
[~, maxIdx] = max(avg_cphyto_values);
[~, minIdx] = min(avg_cphyto_values);
hold on;
plot(timeVector(maxIdx), avg_cphyto_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_cphyto_values(maxIdx), ...
    sprintf('Max: %.2f', avg_cphyto_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_cphyto_values(minIdx), 'ro');
text(timeVector(minIdx), avg_cphyto_values(minIdx), ...
    sprintf('Min: %.2f', avg_cphyto_values(minIdx)), 'VerticalAlignment', 'top');
hold off;

% Save figure
saveas(gcf, 'tit_average_phytobiomass.png');

%% Plotting and Saving PAR
figure;
plot(timeVector, avg_par_values, 'Color', [0.8500 0.4250 0.0980]);
title('tit Average PAR Over Time');
xlabel('Date');
ylabel('Average PAR (mol photons^-^m^2^-^d)');

% Add callouts
[~, maxIdx] = max(avg_par_values);
[~, minIdx] = min(avg_par_values);
hold on;
plot(timeVector(maxIdx), avg_par_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_par_values(maxIdx), ...
    sprintf('Max: %.2f', avg_par_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_par_values(minIdx), 'ro');
text(timeVector(minIdx), avg_par_values(minIdx), ...
    sprintf('Min: %.2f', avg_par_values(minIdx)), 'VerticalAlignment', 'top');
hold off;

% Save figure
saveas(gcf, 'tit_average_par.png');

%% Plotting and Saving SST
figure;
plot(timeVector, avg_sst_values, 'Color', [0.8350 0.0780 0.2840]);
title('Tit Average Sea Surface Temperature Over Time');
xlabel('Date');
ylabel('Average SST (°C)');

% Add callouts
[~, maxIdx] = max(avg_sst_values);
[~, minIdx] = min(avg_sst_values);
hold on;
plot(timeVector(maxIdx), avg_sst_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_sst_values(maxIdx), ...
    sprintf('Max: %.2f°C', avg_sst_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_sst_values(minIdx), 'ro');
text(timeVector(minIdx), avg_sst_values(minIdx), ...
    sprintf('Min: %.2f°C', avg_sst_values(minIdx)), 'VerticalAlignment', 'top');
hold off;

% Save figure
saveas(gcf, 'tit_average_sst.png');

%% Plotting and Saving MLD
figure;
plot(timeVector, avg_mld_values, 'Color', [0.3 0.8 0.3]);
title('Tit Average Mixed Layer Depth Over Time');
xlabel('Date');
ylabel('Average MLD (m)');

% Add callouts
[~, maxIdx] = min(avg_mld_values);
[~, minIdx] = max(avg_mld_values);
hold on;
plot(timeVector(maxIdx), avg_mld_values(maxIdx), 'ro');
text(timeVector(maxIdx), avg_mld_values(maxIdx), ...
    sprintf('Max: %.2f m', avg_mld_values(maxIdx)), 'VerticalAlignment', 'bottom');
plot(timeVector(minIdx), avg_mld_values(minIdx), 'ro');
text(timeVector(minIdx), avg_mld_values(minIdx), ...
    sprintf('Min: %.2f m', avg_mld_values(minIdx)), 'VerticalAlignment', 'top');
hold off;

% Save figure
saveas(gcf, 'tit_average_mld.png');

%% Create a table with the time and mean values for each variable
T = table(timeVector', avg_chla_values, avg_par_values, avg_sst_values, ...
    avg_bbp443_values, avg_cphyto_values, avg_mld_values, ...
    'VariableNames', {'Time', 'Average_Chlorophyll', 'Average_PAR', 'Average_SST', ...
    'Average_BBP443', 'Average_Cphyto', 'Average_MLD'});

% Save the table as an Excel file
writetable(T, 'tit_average_values_over_time.xlsx');