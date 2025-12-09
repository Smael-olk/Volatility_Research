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


%% Hedging Strategy 
% Automates the Hedging Strategy Backtest (Question 6)
% Dates: 08 Dec (Start), 11 Dec, 12 Dec, 13 Dec. 

% --- 1. CONFIGURATION ---
Dates = {'2017-12-08', '2017-12-11', '2017-12-12', '2017-12-13'};
Results_File = './params/CalibratedParams.mat';
FinalDate = datetime('09-Dec-2019'); % Product Maturity

% Positions
Pos_Cert = -1;   % Short the Certificate
Pos_Underlying = 0; 
Pos_Put_Option = 0; 
Cash_Account   = 0; 

% Storage
Prev_Port_Value = 0;
PnL_History = [];

% Load Data
if ~isfile(Results_File), error('File not found'); end
load(Results_File, 'Results');

fprintf('Starting Strategy: DELTA-VEGA NEUTRAL (Consistent FD Greeks)...\n');

for t = 1:length(Dates)
    CurrentDateStr = Dates{t};
    CurrentDate = datetime(CurrentDateStr, 'InputFormat', 'yyyy-MM-dd');
    
    % --- 2. GET MARKET DATA ---
    found_idx = -1;
    for j = 1:length(Results)
        if isequal(dateshift(Results(j).Date, 'start', 'day'), ...
                   dateshift(CurrentDate, 'start', 'day'))
            found_idx = j; break;
        end
    end
    if found_idx == -1, continue; end
    
    DayData = Results(found_idx);
    params  = DayData.Params; 
    S0      = DayData.Spot;
    
    % Interpolate Rates
    T_final = years(FinalDate - CurrentDate);
    r_curve = -log(DayData.Discounts) ./ DayData.T_years;
    r_riskfree = interp1(DayData.T_years, r_curve, T_final, 'linear', 'extrap');
    F_2y = interp1(DayData.T_years, DayData.Fwds, T_final, 'linear', 'extrap');
    q_div = r_riskfree - log(F_2y / S0) / T_final;
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('DATE: %s | S0: %.2f | r: %.2f%% \n', CurrentDateStr, S0, r_riskfree*100);
    fprintf('--------------------------------------------------\n');
    
    % --- 3. VALUATION (Mark-to-Market) ---
    [Price_Cert, ~] = Price_Miami_Certificate_Function(params, S0, r_riskfree, q_div, CurrentDate);
    Val_Cert = Pos_Cert * Price_Cert;
    
    % Price Existing Hedges
    Val_Option = 0;
    if Pos_Put_Option ~= 0
        T_rem = years(Current_T - CurrentDate);
        if T_rem > 0.005
            Price_Put_Curr = NIG_Pricer(params, S0, Current_K, exp(-r_riskfree*T_rem), T_rem, 2); 
            Val_Option = Pos_Put_Option * Price_Put_Curr;
        else
            Payoff = max(Current_K - S0, 0);
            Cash_Account = Cash_Account + (Pos_Put_Option * Payoff);
            Pos_Put_Option = 0;
        end
    end
    Val_Stock = Pos_Underlying * S0;
    
    % --- 4. CALCULATE P&L ---
    if t > 1
        dt_days = days(CurrentDate - datetime(Dates{t-1}, 'InputFormat', 'yyyy-MM-dd'));
        Interest = Cash_Account * (r_riskfree * dt_days/365);
        Cash_Account = Cash_Account + Interest;
        
        Port_Value = Val_Cert + Val_Stock + Val_Option + Cash_Account;
        Daily_PnL = Port_Value - Prev_Port_Value;
        PnL_History(end+1) = Daily_PnL;
        
        fprintf('   Daily P&L:         $ %+.2f\n', Daily_PnL);
    else
        Cash_Account = -Val_Cert; 
        Prev_Port_Value = 0; 
    end
    Prev_Port_Value = Val_Cert + Val_Stock + Val_Option + Cash_Account;
    
    % --- 5. CALCULATE GREEKS (CONSISTENT METHOD) ---
    dS = S0 * 0.005; dVol = params(1) * 0.005;
    
    % A. Certificate Greeks
    P_up = Price_Miami_Certificate_Function(params, S0+dS, r_riskfree, q_div, CurrentDate);
    P_dn = Price_Miami_Certificate_Function(params, S0-dS, r_riskfree, q_div, CurrentDate);
    Delta_Cert = (P_up - P_dn) / (2*dS);
    
    p_up = params; p_up(1) = params(1) + dVol;
    p_dn = params; p_dn(1) = params(1) - dVol;
    P_v_up = Price_Miami_Certificate_Function(p_up, S0, r_riskfree, q_div, CurrentDate);
    P_v_dn = Price_Miami_Certificate_Function(p_dn, S0, r_riskfree, q_div, CurrentDate);
    Vega_Cert = (P_v_up - P_v_dn) / (2*dVol);
    
    Total_Vega_Risk  = Pos_Cert * Vega_Cert;   % Your SHORT position risk
    Total_Delta_Risk = Pos_Cert * Delta_Cert;
    
    fprintf('   [Risk] Cert (Short): Delta %+.0f | Vega %+.0f\n', Total_Delta_Risk, Total_Vega_Risk);
    
    % --- 6. HEDGE INSTRUMENT GREEKS ---
    T_new = 0.25; K_new = S0; Target_Exp = CurrentDate + days(90);
    
    % Put Greeks
    Put_up = NIG_Pricer(params, S0+dS, K_new, exp(-r_riskfree*T_new), T_new, 2);
    Put_dn = NIG_Pricer(params, S0-dS, K_new, exp(-r_riskfree*T_new), T_new, 2);
    Delta_Put = (Put_up - Put_dn) / (2*dS);
    
    Put_v_up = NIG_Pricer(p_up, S0, K_new, exp(-r_riskfree*T_new), T_new, 2);
    Put_v_dn = NIG_Pricer(p_dn, S0, K_new, exp(-r_riskfree*T_new), T_new, 2);
    Vega_Put = (Put_v_up - Put_v_dn) / (2*dVol);
    
    % --- 7. SOLVE SYSTEM ---
    Target_Pos_Put = -Total_Vega_Risk / Vega_Put;
    Target_Pos_Stock = -(Total_Delta_Risk + Target_Pos_Put * Delta_Put);
    
    % --- 8. VERIFICATION & DISPLAY ---
    % Calculate Risk of the Hedges
    Hedge_Delta = (Target_Pos_Stock * 1) + (Target_Pos_Put * Delta_Put);
    Hedge_Vega  = (Target_Pos_Put * Vega_Put);
    
    % Calculate Net Position (Should be ~0)
    Net_Delta = Total_Delta_Risk + Hedge_Delta;
    Net_Vega  = Total_Vega_Risk + Hedge_Vega;
    
    fprintf('   [Hedge] Allocation:  %.0f Puts    | %.0f Stock\n', Target_Pos_Put, Target_Pos_Stock);
    fprintf('   [Net]   Post-Hedge:  Delta %+.2f  | Vega %+.2f  <-- (Target: 0)\n', Net_Delta, Net_Vega);
    
    % --- 9. EXECUTE ---
    if Pos_Put_Option ~= 0
        Cost_Close = -Pos_Put_Option * Price_Put_Curr; 
        Cash_Account = Cash_Account - Cost_Close;
    end
    
    Price_New_Put = NIG_Pricer(params, S0, K_new, exp(-r_riskfree*T_new), T_new, 2);
    Cost_Open = Target_Pos_Put * Price_New_Put;
    Cash_Account = Cash_Account - Cost_Open;
    
    Trade_Stock = Target_Pos_Stock - Pos_Underlying;
    Cash_Account = Cash_Account - (Trade_Stock * S0);
    
    Pos_Put_Option = Target_Pos_Put;
    Pos_Underlying = Target_Pos_Stock;
    Current_K = K_new;
    Current_T = Target_Exp;
end
fprintf('\nFinal P&L (Mean Squared Error): %.4f\n', mean(PnL_History.^2));

%% --- PLOTTING P&L EVOLUTION (BLOOMBERG STYLE) ---

% 1. Prepare Data
PlotDates = datetime(Dates, 'InputFormat', 'yyyy-MM-dd');

if length(PnL_History) == length(PlotDates) - 1
    PnL_Aligned = [0, PnL_History];
    SqPnL_Aligned = [0, PnL_History.^2];
else
    PnL_Aligned = PnL_History;
    SqPnL_Aligned = PnL_History.^2;
end
CumPnL = cumsum(PnL_Aligned);

% --- STYLE DEFINITIONS ---
bbg_black  = [0 0 0];          % Background
bbg_orange = [1.0 0.6 0.0];    % "Amber" Data Color
bbg_white  = [1 1 1];          % Text Color
bbg_grid   = [0.3 0.3 0.3];    % Subtle Grid Lines

% 2. Create Figure
figure('Name', 'Strategy Performance Analysis', ...
       'Color', bbg_black, ... % Black Figure Background
       'Position', [100 100 1000 800], ...
       'InvertHardcopy', 'off');

% --- Subplot 1: Cumulative P&L ---
ax1 = subplot(2,1,1);
set(ax1, 'Color', bbg_black, 'XColor', bbg_white, 'YColor', bbg_white, ...
    'GridColor', bbg_grid, 'GridAlpha', 0.6); % Axes Styling

hold on;
plot(PlotDates, CumPnL, '-o', ...
    'Color', bbg_orange, ...
    'LineWidth', 2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', bbg_orange, ...
    'MarkerEdgeColor', bbg_orange);

yline(0, '--', 'Color', bbg_white, 'LineWidth', 1.5); % White break-even line

title('Cumulative P&L Evolution (Net Wealth)', 'Color', bbg_white);
ylabel('USD', 'Color', bbg_white);
grid on;
grid minor;

% Annotate Final Value (White Text)
text(PlotDates(end), CumPnL(end), sprintf('  $%.0f', CumPnL(end)), ...
    'VerticalAlignment', 'bottom', ...
    'Color', bbg_white, ... 
    'FontSize', 10, 'FontWeight', 'bold');

% --- Subplot 2: Squared P&L ---
ax2 = subplot(2,1,2);
set(ax2, 'Color', bbg_black, 'XColor', bbg_white, 'YColor', bbg_white, ...
    'GridColor', bbg_grid, 'GridAlpha', 0.6); % Axes Styling

hold on;
b = bar(PlotDates, SqPnL_Aligned, ...
    'FaceColor', bbg_orange, ...
    'EdgeColor', 'none'); 

title('Daily P&L Squared (Contribution to Error Metric)', 'Color', bbg_white);
ylabel('USD ^2', 'Color', bbg_white);
grid on;

% Highlight the worst day
[max_sq_err, idx_max] = max(SqPnL_Aligned);
text(PlotDates(idx_max), max_sq_err, '  Worst Day', ...
    'VerticalAlignment', 'bottom', ...
    'Color', bbg_white, ... % White text for annotation
    'FontWeight', 'bold');

% --- Global Title ---
Metric_Plot = mean(PnL_History.^2); 
sgt = sgtitle(['Strategy Results | RMSE: $' num2str(sqrt(Metric_Plot), '%.2f')], ...
    'Color', bbg_white, ...
    'FontSize', 14, 'FontWeight', 'bold');