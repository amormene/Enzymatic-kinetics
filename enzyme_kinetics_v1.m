%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ENZYME_KINETICS_GUI.m
%
% AUTHOR: Andrés Mormeneo-Segarra
% MAIL: amormene@uji.es; anmorse@upv.es
% DATE: June 2025
%
% DESCRIPTION:
% -------------------------------------------------------------------------
% This MATLAB GUI allows the user to load experimental data from an 
% enzyme kinetics experiment (initial substrate concentration vs. 
% initial reaction rate), and automatically:
%   - Asks for concentration and time units
%   - Computes the Michaelis-Menten parameters (K_M and r_m)
%     using four linearization models:
%        1. Michaelis-Menten (nonlinear regression)
%        2. Lineweaver-Burk
%        3. Eadie-Hofstee
%        4. Hanes-Woolf
%   - Displays a figure with all four models and fitted curves
%   - Displays a table of kinetic parameters with proper units
%   - Allows saving of processed data with all linearized transformations
%
% INPUT FORMAT:
% -------------------------------------------------------------------------
% A plain text file (*.txt) with two columns:
%   Column 1: Initial substrate concentrations [c_S0]
%   Column 2: Initial rates [r_0]
% No headers should be included.
%
% USER INTERACTION:
% -------------------------------------------------------------------------
% 1. Press "Load Data" to select your file.
% 2. Enter:
%    - Units of concentration (e.g., 'mol/L', 'mM')
%    - Units of time (e.g., 's', 'min', 'h')
%    * Velocity units are automatically assumed as [concentration/time]
% 3. All figures and results are automatically shown.
% 4. A processed data table can be exported as a .txt file.
%
% REQUIREMENTS:
% -------------------------------------------------------------------------
% MATLAB version: R2014b or later recommended
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function enzyme_kinetics_gui()
    clc; clear all; clearvars; close all;
    f = figure('Name', 'Enzyme Kinetics Problem', 'NumberTitle', 'off', ...
               'Position', [500, 300, 400, 150]);
    uicontrol('Style', 'pushbutton', 'String', 'Load Data', ...
              'Position', [100 50 200 40], 'FontSize', 12, ...
              'Callback', @load_data);
end

function load_data(~, ~)
    [file, path] = uigetfile('*.txt', 'Select Data File');
    if isequal(file, 0)
        return;
    end
    data = readmatrix(fullfile(path, file));
    c_S0 = data(:, 1);
    r_0 = data(:, 2);

    prompt = {'Concentration units (e.g., mol/L):', 'Time units (e.g., s, min, h):'};
    dlg_title = 'Input Units';
    dims = [1 35];
    definput = {'mol/L', 's'};
    answer = inputdlg(prompt, dlg_title, dims, definput);

    if isempty(answer)
        return;
    end

    conc_unit = answer{1};
    time_unit = answer{2};
    rate_unit = [conc_unit '/' time_unit];

    assignin('base', 'conc_unit', conc_unit);
    assignin('base', 'time_unit', time_unit);
    assignin('base', 'rate_unit', rate_unit);

    plot_initial_data(c_S0, r_0, conc_unit, rate_unit);
    linearize_data(c_S0, r_0, conc_unit, rate_unit, time_unit);
    save_processed_data(c_S0, r_0);
end

function plot_initial_data(c_S0, r_0, conc_unit, rate_unit)
    figure('Name','Initial Rate vs Concentration');
    scatter(c_S0, r_0, 'filled', 'MarkerFaceColor', 'b');
    xlabel(['c_{S0} (' conc_unit ')']);
    ylabel(['r_0 (' rate_unit ')']);
    title('Initial Rate vs Concentration');
    grid on;
end

function linearize_data(c_S0, r_0, conc_unit, rate_unit, time_unit)
    x_lb = 1 ./ c_S0;
    y_lb = 1 ./ r_0;
    [slb, ilb, r2_lb] = fit_linear(x_lb, y_lb);
    Vmax_lb = 1 / ilb;
    Km_lb = slb * Vmax_lb;

    x_eh = r_0 ./ c_S0;
    y_eh = r_0;
    [seh, ieh, r2_eh] = fit_linear(x_eh, y_eh);
    Vmax_eh = ieh;
    Km_eh = -seh;

    x_hw = c_S0;
    y_hw = c_S0 ./ r_0;
    [shw, ihw, r2_hw] = fit_linear(x_hw, y_hw);
    Vmax_hw = 1 / shw;
    Km_hw = ihw * Vmax_hw;

    mm_model = @(b, x)(b(1) * x ./ (b(2) + x));
    beta0 = [max(r_0), mean(c_S0)];
    opts = statset('nlinfit');
    beta = nlinfit(c_S0, r_0, mm_model, beta0, opts);
    Vmax_nl = beta(1);
    Km_nl = beta(2);

    plot_all(c_S0, r_0, ...
             x_lb, y_lb, slb, ilb, ...
             x_eh, y_eh, seh, ieh, ...
             x_hw, y_hw, shw, ihw, ...
             Vmax_nl, Km_nl, conc_unit, rate_unit, time_unit);

    show_results(slb, ilb, r2_lb, Vmax_lb, Km_lb, ...
                 seh, ieh, r2_eh, Vmax_eh, Km_eh, ...
                 shw, ihw, r2_hw, Vmax_hw, Km_hw, ...
                 Vmax_nl, Km_nl, conc_unit, time_unit);
end

function [slope, intercept, r2] = fit_linear(x, y)
    X = [x(:) ones(length(x),1)];
    coeffs = X \ y(:);
    slope = coeffs(1);
    intercept = coeffs(2);
    y_fit = X * coeffs;
    SS_res = sum((y - y_fit).^2);
    SS_tot = sum((y - mean(y)).^2);
    r2 = 1 - SS_res / SS_tot;
end

function save_processed_data(c_S0, r_0)
    output = [c_S0, r_0, 1./r_0, 1./c_S0, r_0./c_S0, c_S0./r_0];
    headers = {'c_S0', 'r_0', '1/r_0', '1/c_S0', 'r_0/c_S0', 'c_S0/r_0'};
    
    [file, path, idx] = uiputfile({'*.txt','Text File (*.txt)'; ...
                                   '*.csv','CSV File (*.csv)'; ...
                                   '*.xlsx','Excel File (*.xlsx)'}, ...
                                   'Save Processed Data');
    
    if isequal(file, 0)
        return;
    end
    
    outpath = fullfile(path, file);

    switch idx
        case 1  % TXT
            fid = fopen(outpath, 'w');
            fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', headers{:});
            fclose(fid);
            dlmwrite(outpath, output, '-append', 'delimiter', '\t', 'precision', '%.6f');
        case 2  % CSV
            data_cell = [headers; num2cell(output)];
            writecell(data_cell, outpath);
        case 3  % XLSX
            data_cell = [headers; num2cell(output)];
            writecell(data_cell, outpath);
    end
end


function plot_all(c_S0, r_0, x_lb, y_lb, slb, ilb, ...
                  x_eh, y_eh, seh, ieh, ...
                  x_hw, y_hw, shw, ihw, ...
                  Vmax_nl, Km_nl, conc_unit, rate_unit, time_unit)

    figure('Name','All Kinetic Plots', 'WindowState', 'maximized');

    % 1. Michaelis-Menten
    subplot(2,2,1);
    scatter(c_S0, r_0, 'b'); hold on;
    xfit = linspace(min(c_S0), max(c_S0), 100);
    yfit = Vmax_nl * xfit ./ (Km_nl + xfit);
    plot(xfit, yfit, 'r', 'LineWidth', 2);
    title('Michaelis-Menten');
    xlabel(['c_{S0} (' conc_unit ')']);
    ylabel(['r_0 (' rate_unit ')']);
    grid on;
    eq1 = sprintf('r_0 = %.3f·c_{S0} / (%.3f + c_{S0})', Vmax_nl, Km_nl);
    text(0.95, 0.10, eq1, 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'Units', 'normalized');

    % 2. Lineweaver-Burk
    subplot(2,2,2);
    scatter(x_lb, y_lb, 'b'); hold on;
    plot(x_lb, slb*x_lb + ilb, 'r', 'LineWidth', 2);
    title('Lineweaver-Burk');
    xlabel(['1/c_{S0} (' conc_unit '^{-1})']);
    ylabel(['1/r_0 (' rate_unit '^{-1})']);
    grid on;
    eq2 = sprintf('1/r_0 = %.3f·(1/c_{S0}) + %.3f', slb, ilb);
    text(0.95, 0.10, eq2, 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'Units', 'normalized');

    % 3. Eadie-Hofstee
    subplot(2,2,3);
    scatter(x_eh, y_eh, 'b'); hold on;
    plot(x_eh, seh*x_eh + ieh, 'r', 'LineWidth', 2);
    title('Eadie-Hofstee');
    xlabel(['r_0 / c_{S0} (' time_unit '^{-1})']);
    ylabel(['r_0 (' rate_unit ')']);
    grid on;
    eq3 = sprintf('r_0 = %.3f·(r_0/c_{S0}) + %.3f', seh, ieh);
    text(0.05, 0.10, eq3, 'FontSize', 10, ...
         'HorizontalAlignment', 'left', 'Units', 'normalized');

    % 4. Hanes-Woolf
    subplot(2,2,4);
    scatter(x_hw, y_hw, 'b'); hold on;
    plot(x_hw, shw*x_hw + ihw, 'r', 'LineWidth', 2);
    title('Hanes-Woolf');
    xlabel(['c_{S0} (' conc_unit ')']);
    ylabel(['c_{S0}/r_0 (' time_unit ')']);
    grid on;
    eq4 = sprintf('c_{S0}/r_0 = %.3f·c_{S0} + %.3f', shw, ihw);
    text(0.95, 0.10, eq4, 'FontSize', 10, ...
         'HorizontalAlignment', 'right', 'Units', 'normalized');
end

function show_results(lb_slope, lb_intercept, r2_lb, Vmax_lb, Km_lb, ...
                      eh_slope, eh_intercept, r2_eh, Vmax_eh, Km_eh, ...
                      hw_slope, hw_intercept, r2_hw, Vmax_hw, Km_hw, ...
                      Vmax_nl, Km_nl, conc_unit, time_unit)

    rate_unit = [conc_unit '/' time_unit];

    data = {
        'Lineweaver-Burk',     sprintf('%.3f', lb_slope),   sprintf('%.3f', lb_intercept), sprintf('%.3f', r2_lb),   sprintf('%.3f', Vmax_lb),   sprintf('%.3f', Km_lb);
        'Eadie-Hofstee',     sprintf('%.3f', eh_slope),   sprintf('%.3f', eh_intercept), sprintf('%.3f', r2_eh),   sprintf('%.3f', Vmax_eh),   sprintf('%.3f', Km_eh);
        'Hanes-Woolf',     sprintf('%.3f', hw_slope),   sprintf('%.3f', hw_intercept), sprintf('%.3f', r2_hw),   sprintf('%.3f', Vmax_hw),   sprintf('%.3f', Km_hw);
        'Average', '-',                         '-',                           '-',                       sprintf('%.3f', mean([Vmax_lb, Vmax_eh, Vmax_hw])), ...
                                                                                                             sprintf('%.3f', mean([Km_lb, Km_eh, Km_hw]));
        'Nonlinear Regression',     '-',                         '-',                           '-',                       sprintf('%.3f', Vmax_nl),   sprintf('%.3f', Km_nl)
    };

    colNames = {
        'Model', ...
        'Slope', ...
        'Intercept', ...
        'R²', ...
        ['rₘ (' rate_unit ')'], ...
        ['Kₘ (' conc_unit ')']};

    fig = figure('Name', 'Kinetic Parameters Table', ...
                 'NumberTitle', 'off', ...
                 'Color', [1 1 1], ...
                 'Position', [300 300 750 220]);

    t = uitable(fig, ...
        'Data', data, ...
        'ColumnName', colNames, ...
        'ColumnEditable', false(1,6), ...
        'ColumnWidth', {200, 90, 90, 60, 130, 130}, ...
        'FontSize', 12, ...
        'FontName', 'Arial', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.1 0.96 0.85]);
end
