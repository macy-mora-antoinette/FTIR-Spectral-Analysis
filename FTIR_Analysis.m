%% FTIR Analysis

clear all
close all
clc

%%
%{
Updated 13 September 2023
Author: Macy Mora-Antoinette

MATLAB R2023a  

This matlab code is intended to help users compute the compositional
properties from FTIR Absorption spectra in bones. 

Reference:
Hunt, H. B. et al. Altered Tissue Composition, Microarchitecture, and 
Mechanical Performance in Cancellous Bone From Men With Type 2 Diabetes 
Mellitus. Journal of Bone and Mineral Research 34, 1191–1206 (2019).
https://doi.org/10.1002/jbmr.3711

It is separated into three parts:

PART 1: A simple for loop to iterate through raw data files (csv)
in a single directory and take the average of all replicates. 
Plots are produced to visualize the average. 

PART 2: Processing the data using the average. Computes the area and peak
values for the different compositional properties. It will save plots for each
area and peak evaluated. Plots are to help the user evaluate data quality.  

PART 3: These functions used in PART 2 to help compute the different
values: spectra_peaks, spectra_area, and baseline.

Key things to keep in mind:
1. Works best if you run in sections PART 1 followed by PART 2
2. This code requires a specific naming format for your csv files. It
expects the sample replicates to have identical names except for #.csv at
the end to indicate the replicate.
Ex filename1.csv, filename2.csv, filename3.csv
2. It does not deal with corner cases, so you may need to manually change
wavelengths inputs for areas/peaks as needed
3. If you want to add additional area or peak computations, make sure you
change the "outputs" variable on line 134 and add the necessary columns to
the table on line 143

                                --- PART 1 ---
%}

%% Get Averages

%select the directory with all the raw data .csv files saved
disp('Select the directory with raw data'); all_files = uigetdir(); 
file_names = dir(all_files);
file_names = {file_names.name};
file_names = string(file_names(3:end))';

%make sure sample replicates have identical names with #.csv at the end 
samples = unique(regexprep(file_names,'\d.csv',''));

%Iterate over all sample names
for i = 1:length(samples)
    sample = samples(i);
    idx = find(contains(file_names,sample));
    
    %make amatrix with all replicates and plot
    avg = [];
    for j = 1:length(idx)
        file = file_names(idx(j));
        temp = readmatrix(append(all_files,'\', file));
        avg(:,j) = temp(:, 2);
        plot(temp(:, 1), temp(:, 2), 'DisplayName',file)
        hold on
    end

    %save data
    dat = [temp(:, 1), mean(avg, 2)];
    T = array2table(dat);
    T.Properties.VariableNames(1:2) = {'inverse_cm','A'};
    writetable(T, append(all_files,'\', sample, ' Average.csv'));

    %plot aesthetics
    plot(dat(:, 1), dat(:, 2), '--k', 'DisplayName', 'Average')
    legend
    title(sample); xlabel('Wavelength (cm^{-1})'); ylabel('Absorbance');
    hold off
    savefig(append(all_files,'\', sample, ' Average.fig'))
end

close all
%%
%{
                                --- PART 2 ---
%% Process Data
%}


%Select the directory where all the averaged data is to be used for analysis
disp('Select the directory with averaged data'); ls_dir = uigetdir();
ls = dir(ls_dir);
ls = {ls.name};
ls = ls(contains(ls, 'Average.csv'))';
samples = regexprep(ls,' Average.csv','');

%Iterate over all average files
Results = {};
for k = 2%1:length(ls)
    file = ls{k};
    sample = samples{k};
    raw = readmatrix(append(ls_dir, '\', file));      

    %Crop to only 800 to 1800 cm-1
    mask = raw(:,1) <=1800;
    df = raw(mask, :);

    %Prep plots
    tiledlayout(2,3)
    
    %Get Baseline
    base  = baseline(df, {[800 850],[1200 1400],[1750 1800]}, 2, sample);
    %base = baseline(df, {[820 840],[880 900], [1300 1400],[1750 1800]}, 3, sample); %use only if you need it
    
    % Get Areas
    phosphArea = spectra_area(base, 916, 1180, sample, 'Phosphate'); %phosphate 890cm-1 to 1200cm-1    
    amideArea = spectra_area(base, 1597, 1712, sample, 'Amide'); %amide 1597cm-1 to 1730cm-1
    carbArea = spectra_area(base, 852, 890, sample, 'Carbonate'); %carbonate 820cm-1 to 890cm-1

    % Get Peaks 
    Peaks_phosphate = spectra_peaks(base, 890, 1200, [1030 1020 1127 1096], sample, 'Phosphate'); %Phosphate Baselined Peaks  
    Peaks_amide = spectra_peaks(base, 1597, 1730, [1660 1690], sample, 'Amide'); %Amide Baselined Peaks

    % Save Results
    outputs = {sample, phosphArea, amideArea, carbArea,...
        Peaks_phosphate(1), Peaks_phosphate(2), Peaks_phosphate(3), Peaks_phosphate(4),...
        Peaks_amide(1), Peaks_amide(2)};
    Results = [Results;outputs];
    savefig(append(ls_dir, '\', sample, ' Processed.fig'))  
end

%store results as a table
T = cell2table(Results);
T.Properties.VariableNames = {'Sample', 'phosphArea', 'amideArea', 'carbArea',...
        'peak1030', 'peak1020', 'peak1127', 'peak1096', 'peak1660', 'peak1690'};
writetable(T, append(ls_dir, '\FTIR_Data.csv'));
close all
%%
%{
                                --- PART 3 ---
%}

%%
function Peaks = spectra_peaks(data, lowerBound,upperBound, peak_wavelengths, sampleName, peakName)    

    % Description: computes the peaks values at given wavelengths
    % Inputs:
    % data - spectral data
    % lowerBound - wavelength for lower bound of spectral curve 
    % upperBound - wavelength for  upper bound of spectral curve
    % peak_wavelenghts - array of wavelengths for peaks
    % sampleName - to name plot
    % peakName - to name plot
    %
    % Outputs:
    % Peaks - array. column 1 is the wavelength, column 2 is the peak value 

    % Baseline
    dat = data( (data(:,1) >= lowerBound) & (data(:,1) <= upperBound), :);
    p = polyfit(dat([1, end], 1), dat([1, end], 2), 1);
    y = polyval(p, dat(:, 1));
    ypoly = dat(:, 2) - y;
    
    %Get Peaks
    Peaks = [];
    nexttile
    hold on
    for i = 1:length(peak_wavelengths)
        wavelength = peak_wavelengths(i);
        peak = ypoly( dat(:,1) == wavelength );
        Peaks(i) = peak;
        plot(wavelength, peak, '*', 'DisplayName', append(num2str(wavelength), ' peak'));
    end
    
    %Plot
    plot(data(:, 1), data(:, 2), 'DisplayName', 'Baselined');
    xline(lowerBound, 'k--', 'DisplayName', append(num2str(lowerBound),'cm^{-1}'));
    xline(upperBound, 'k--', 'DisplayName', append(num2str(upperBound),'cm^{-1}'));
    plot(dat(:, 1), ypoly, 'DisplayName', append( peakName, ' Baseline'));
    %aesthetics
    xlin = linspace(dat(1, 1), dat(end, 1), length(data));
    ylin = linspace(y(1), y(end), length(data));
    plot(xlin, ylin, 'DisplayName', 'Poly n=1');
    title(append(sampleName,' ', peakName, ' Peaks')); xlabel('Wavelength (cm^{-1})'); ylabel('Absorbance');
    legend
    hold off
end

%%
function Area = spectra_area(data, lowerBound,upperBound, sampleName, peakName)

    % Description: computes the areas at given wavelengths
    % Inputs:
    % data - spectral data
    % lowerBound - wavelength for lower bound of spectral curve 
    % upperBound - wavelength for  upper bound of spectral curve
    % sampleName - to name plot
    % peakName - to name plot
    %
    % Outputs:
    % Area - scalar. area value    

    % Baseline 
    dat = data( (data(:,1) >= lowerBound) & (data(:,1) <= upperBound), :);
    p = polyfit(dat([1, end], 1), dat([1, end], 2), 1);
    y = polyval(p, dat(:, 1));
    ypoly = dat(:, 2) - y;
    
    %Compute Area within bounds
    Area = trapz(ypoly);
    
    %Plot
    nexttile
    plot(data(:, 1), data(:, 2), 'DisplayName', 'Baselined');
    hold on
    xline(lowerBound, 'k--', 'DisplayName', append(num2str(lowerBound),'cm^{-1}'));
    xline(upperBound, 'k--', 'DisplayName', append(num2str(upperBound),'cm^{-1}'));
    area(dat(:, 1), ypoly, 'DisplayName', append( peakName, ' Area'));
    %aesthetics
    xlin = linspace(dat(1, 1), dat(end, 1), length(data));
    ylin = linspace(y(1), y(end), length(data));
    plot(xlin, ylin, 'DisplayName', 'Poly n=1');
    title(append(sampleName,' ', peakName, ' Area')); xlabel('Wavelength (cm^{-1})'); ylabel('Absorbance');
    legend
    hold off
end

%%
function dat_baselined = baseline(data, limits, degree, sampleName)
    
    % Description: computes baseline with local minima
    % Inputs:
    % data - spectral data
    % limits - cell array with limits for local minima 
    % degree - degree of polynomial baseline. Typically 2. 3 for rare cases
    % sampleName - to name plot
    %
    % Outputs:
    % dat_baselined    

    %Get minima
    mins = [];
    for i=1:length(limits)
        lims = limits{i};
        val = local_min(data, lims(1), lims(2));
        mins(i, :) = val(end, :);
    end

    %Get polynomial
    p = polyfit(mins(:, 1),mins(:, 2),degree);
    y = polyval(p,data(:,1));
    
    %Plot
    nexttile
    plot(data(:,1), data(:, 2), 'DisplayName', 'Raw');
    hold on
    plot(data(:,1), y, 'DisplayName', append('poly n=', num2str(degree)));
    base = data(:, 2) - y;
    plot(data(:,1), base, 'DisplayName', 'Baselined');
    title(append(sampleName, ' Baselined')); xlabel('Wavelength (cm^{-1})'); ylabel('Absorbance');
    legend
    hold off
    dat_baselined = [data(:,1), base];

    %% find local minima
    function minimum = local_min(data, lowerLim, upperLim)
        range = data( (data(:,1) >= lowerLim) & (data(:,1) <= upperLim), :);
        localMin = min(range);
        minimum = range(range(:, 2) == localMin(2), :);
    end
end