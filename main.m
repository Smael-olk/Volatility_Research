clc
close all
clear
warning('off','all')
format long
rng(42)
path(pathdef)
%% add required folders
addpath("Code\");
h_subfolders = genpath('Code');
addpath(h_subfolders);
clear h_subfolders;
addpath("Dati_Train");
addpath("Plots\");
addpath("Lewis\")
addpath("plotting_functions\")
addpath("params\")


%% load data
filelist = dir('.\Dati_Train\*.csv');
numFiles = length(filelist);
spot_file_name = 'spot_SPX_hist.csv';

%% Question a : Term Structure on 08 Dec 2017
target_date = '2017-12-08';
% Find the specific file for this date
file_match = dir(fullfile('Dati_Train', ['*' target_date '*.csv']));

if isempty(file_match)
    error('File for date %s not found.', target_date);
end
full_path = fullfile(file_match(1).folder, file_match(1).name);
% Call IR_curves with the actual file path
[discounts,forward_prices,expiries,TradeDate,Spot,Option_data] = IR_curves(full_path, spot_file_name,0.1, 0.6, 5, 200);
T_years = years(expiries - TradeDate);
R_years = -log(discounts(:)') ./ T_years(:)';
% Plotting
plot_term_struct(discounts,R_years,forward_prices,TradeDate,Spot,T_years);

%% Running calibration on whole 6 months period : Already Done but can be repeated for Verification.
run_calibration();
%% Question b : Volatility Surface 08-Dec-2017
% 1. Configuration
target_date_str = '2017-12-08'; % The specific date you want to check
results_file    = './params/CalibratedParams.mat';

% 2. Load the Big Results File
if ~isfile(results_file)
    error('Results file not found: %s. \nPlease run run_calibration.m first.', results_file);
end

fprintf('Loading results from %s...\n', results_file);
load(results_file, 'Results');

% 3. Search for the Target Date
target_date_dt = datetime(target_date_str, 'InputFormat', 'yyyy-MM-dd');
found_idx = -1;

for j = 1:length(Results)
    if isequal(dateshift(Results(j).Date, 'start', 'day'), dateshift(target_date_dt, 'start', 'day'))
        found_idx = j;
        break;
    end
end

if found_idx == -1
    error('Date %s not found in CalibratedParams.mat.', target_date_str);
end

fprintf('Found data for %s at index %d.\n', target_date_str, found_idx);

% 4. Extract Data for that Day
DayData = Results(found_idx);

if isempty(DayData.Params)
    error('Calibration parameters are empty for this date.');
end

optimal_params = DayData.Params;
Spot           = DayData.Spot;
TradeDate      = DayData.Date;
T_years        = DayData.T_years;
discounts      = DayData.Discounts;
fwds           = DayData.Fwds;

fprintf('Params: %s\n', mat2str(optimal_params, 3));

% 5. Regenerate Dense Grid
% We populate specific columns (1,2,3,5) so arbitrage_check can read them.
num_expiries = length(T_years);
plot_data = cell(num_expiries, 8); 

tol = 1e-6; 
max_iter = 200;

fprintf('Regenerating surface data...\n');

for i = 1:num_expiries
    Dt = T_years(i);
    B  = discounts(i);
    F0 = fwds(i);
    
    if Dt <= 0, continue; end
    
    % A. Create Dense Strike Grid (Moneyness 80% to 120%)
    K_dense = linspace(Spot * 0.80, Spot * 1.20, 50)';
    
    % B. Calculate NIG Model Prices (Calls)
    C_model = NIG_Pricer(optimal_params, F0, K_dense, B, Dt);
    
    % C. Invert to get Implied Volatility
    iv_model = zeros(size(C_model));
    for k = 1:length(C_model)
        try
            iv_model(k) = Black76Inverse(C_model(k), F0, K_dense(k), B, Dt, tol, max_iter);
        catch
            iv_model(k) = NaN;
        end
    end
    
    % D. Store in Cell Array
    % KEY CHANGE: Storing B and F0 so arbitrage_check can use them
    plot_data{i, 1} = Dt;       % Time
    plot_data{i, 2} = B;        % Discount (Required for Arb Check)
    plot_data{i, 3} = F0;       % Forward  (Required for Arb Check)
    plot_data{i, 5} = K_dense;  % Strikes  (Required for Arb Check limits)
    
    plot_data{i, 7} = iv_model; % IV       (Required for Plot)
    plot_data{i, 8} = K_dense;  % Strikes  (Required for Plot)
end

% 6. Plot Volatility Surface
fprintf('Plotting Surface...\n');
fig = figure('Name', ['NIG Surface: ' datestr(TradeDate)], 'Color', 'k'); 
main_ax = axes('Parent', fig); 
vol_surface_plot(plot_data, main_ax, TradeDate, Spot);

% 7. Run Arbitrage Check
fprintf('Running Arbitrage Check...\n');
arbitrage_check(plot_data, optimal_params);
disp('Analysis Complete.');

% 8. Display Calibration Error Metrics
% Since you have already run the calibration, these values are stored in DayData.
if isfield(DayData, 'RMSE') && ~isempty(DayData.RMSE)
    fprintf('\n================================================\n');
    fprintf('CALIBRATION ACCURACY (%s)\n', datestr(TradeDate));
    fprintf('================================================\n');
    fprintf('RMSE (Root Mean Square Error):       %.4f\n', DayData.RMSE);
    fprintf('MAPE (Mean Absolute Percent Error):  %.2f%%\n', DayData.MAPE * 100);
    fprintf('------------------------------------------------\n');
else
    fprintf('\n[Warning] RMSE/MAPE not found in Results struct.\n');
    fprintf('Ensure you ran the updated run_calibration.m that saves these metrics.\n');
end


%% Question c : Monte Carlo Simulation 
MCSimulation();
compare_pricing(DayData, 5, 1e6, 14);

%% PRICE_MIAMI_CERTIFICATE.m
% 1. Configuration
target_date_str = '2017-12-08'; % The specific date you want to price
results_file    = './params/CalibratedParams.mat';

% 2. Load the Big Results File
if ~isfile(results_file)
    error('Results file not found: %s. \nPlease run run_calibration.m first.', results_file);
end

fprintf('Loading results from %s...\n', results_file);
load(results_file, 'Results');

% 3. Search for the Target Date
target_date_dt = datetime(target_date_str, 'InputFormat', 'yyyy-MM-dd');
found_idx = -1;

for j = 1:length(Results)
    if isequal(dateshift(Results(j).Date, 'start', 'day'), dateshift(target_date_dt, 'start', 'day'))
        found_idx = j;
        break;
    end
end

if found_idx == -1
    error('Date %s not found in CalibratedParams.mat.', target_date_str);
end

fprintf('Found data for %s at index %d.\n', target_date_str, found_idx);

% 4. Extract Data for that Day
DayData = Results(found_idx);

if isempty(DayData.Params)
    error('Calibration parameters are empty for this date.');
end

optimal_params = DayData.Params;
Spot           = DayData.Spot;
TradeDate      = DayData.Date;
T_years        = DayData.T_years;
discounts      = DayData.Discounts;
fwds           = DayData.Fwds; % Corrected: Extracted as 'fwds'

fprintf('Params: %s\n', mat2str(optimal_params, 3));

% 5. Certificate Setup
S0 = Spot;
FinalDate = datetime('09-Dec-2019');
ValuationDate = TradeDate; % Corrected: Use the loaded TradeDate, not hardcoded
Notional = 15e6; 

% --- 6. Interpolate Rates (r, q) ---
T_final = years(FinalDate - ValuationDate);

% Calculate continuous rates from Discount Factors
% r = -ln(B) / T
r_curve = -log(discounts) ./ T_years;

% Interpolate risk-free rate for the specific maturity
r_riskfree = interp1(T_years, r_curve, T_final, 'spline', 'extrap');

% Interpolate Forward Price (Use 'fwds', not 'fwd')
F_2y = interp1(T_years, fwds, T_final, 'spline', 'extrap');

% Calculate Dividend Yield (q) from Forward Parity
% F = S * exp((r-q)T)  =>  q = r - ln(F/S)/T
q_div = r_riskfree - log(F_2y / S0) / T_final;

fprintf('Pricing Parameters:\n r = %.4f%%\n q = %.4f%%\n', r_riskfree*100, q_div*100);

% --- 7. Run Monte Carlo Pricing ---
% Ensure you have the function 'Price_Miami_Certificate_Function' in your path
fprintf('Running Monte Carlo Simulation...\n');

[Price_Cert, StdErr] = Price_Miami_Certificate_Function(optimal_params, S0, r_riskfree, q_div, ValuationDate);

% --- 8. Display Results ---
fprintf('\n--------------------------------------\n');
fprintf('MIAMI CROCODILE CERTIFICATE PRICE (DISCRETE BARRIER)\n');
fprintf('--------------------------------------\n');
fprintf('Price (Total):   $ %.2f\n', Price_Cert); % Assuming function returns total value
fprintf('%% of Notional:   %.2f %%\n', (Price_Cert/Notional)*100);
fprintf('Std Error:       $ %.2f\n', StdErr);
fprintf('--------------------------------------\n');

%% Hedging strategy 6 months from 8 june to 08 december 2017.
% Automates the Hedging Strategy Backtest (6-Month Stress Test)
% Strategy: Short Risk Reversal (Delta-Vega-Gamma Neutral)
% Period: 08 June 2017 to 08 Dec 2017

% --- 1. CONFIGURATION ---
Results_File = './params/CalibratedParams.mat';
FinalDate = datetime('09-Dec-2019'); % Product Maturity
StartDate = datetime('08-Jun-2017');
EndDate   = datetime('08-Dec-2017');
Hedging(StartDate,FinalDate)
%% --- PERFORMANCE EVALUATION ---

% 1. Calculate Metrics
% Portfolio Value (Cumulative P&L over time)
Portfolio_Value = cumsum(PnL_History); 

% Daily P&L Squared
PnL_Squared = PnL_History.^2;

% PERFORMANCE METRIC: Mean Squared Error (1/T * Sum(PnL^2))
% This is the specific metric requested to determine the "winner"
% Summary Stats
MSE = mean(PnL_History.^2);
Cum_PnL = sum(PnL_History);
fprintf('\n--------------------------------------------------\n');
fprintf('BACKTEST COMPLETE\n');
fprintf('Cumulative P&L: $%.2f\n', Cum_PnL);
fprintf('Mean Squared Error (MSE): %.4f\n', MSE);
% Optional: Plot P&L
bbg_orange = [1.0 0.6 0.0];
figure; plot(Dates_History, PnL_History,'color',bbg_orange); 
title('Daily Hedging P&L (June - Dec)'); xlabel('Date'); ylabel('P&L ($)');
grid on;

%% Hedging 11 december to 13 december 2017
% --- 1. CONFIGURATION ---

FinalDate = datetime('09-Dec-2019'); % Product Maturity
StartDate = datetime('08-Dec-2017');
EndDate   = datetime('13-Dec-2017');
Hedging(StartDate,EndDate,FinalDate);
%%
% Summary Stats
MSE = mean(PnL_History.^2);
Cum_PnL = sum(PnL_History);
fprintf('\n--------------------------------------------------\n');
fprintf('BACKTEST COMPLETE\n');
fprintf('Cumulative P&L: $%.2f\n', Cum_PnL);
fprintf('Mean Squared Error (MSE): %.4f\n', MSE);
% Optional: Plot P&L
bbg_orange = [1.0 0.6 0.0];
figure; plot(Dates_History, PnL_History,'color',bbg_orange); 
title('Daily Hedging P&L (June - Dec)'); xlabel('Date'); ylabel('P&L ($)');
grid on;