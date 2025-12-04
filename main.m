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
[discounts,forward_prices,expiries,TradeDate,Spot,Option_data] = IR_curves(full_path, spot_file_name);T_years = years(expiries - TradeDate);
R_years = -log(discounts(:)') ./ T_years(:)';
% Plotting
plot_term_struct(discounts,R_years,forward_prices,TradeDate,Spot,T_years)

%% Useful function to visualize Term Structure Evolution for further inspection / refining of hedging strategy.
%term_structure_evo;

%% Market Data 

marketdata_08 = cell(length(Option_data.datesExpiry), 8);
for i=1:length(Option_data.datesExpiry)
    marketdata_08{i,1} = T_years(i);
    marketdata_08{i,2} = discounts(i);
    marketdata_08{i,3} = forward_prices(i);
    marketdata_08{i,4} = (Option_data.callAsk(i).prices + Option_data.callBid(i).prices ) / 2; % observed market prices
    marketdata_08{i,5} = Option_data.strikes(i).value;
end

%%
sigma = 0.3; % Initial guess (often 0.2 to 0.4 works well)
tol = 1e-6;
max_iter = 100;
for i=1:length(Option_data.datesExpiry)
    iv= [];
    for j=1:length(marketdata_08{i,4})
        T = marketdata_08{i,1};
        B= marketdata_08{i,2};
        F= marketdata_08{i,3};
        C_target = marketdata_08{i,4}(j);
        K = marketdata_08{i,5}(j);
        iv(j)= Black76Inverse(C_target, F, K, B, T, tol, max_iter);
    end
    marketdata_08{i,6} = iv;
end
fig = figure('Name', ' Black76 Market Implied Volatility Surface', 'Position', [100 100 1200 600]);
main_ax = axes('Parent', fig); % Axes handle for plotting
vol_surface_black(marketdata_08,main_ax,TradeDate,Spot);
%% Question b : Volatility Surface 08-Dec-2017


% 1. Configuration & Path Setup
addpath(genpath('Code'));       % Ensure your functions are visible
addpath('Dati_Train');          % Ensure data folder is visible

data_folder     = '.\Dati_Train';
spot_file_name  = 'spot_SPX_hist.csv'; 
target_date     = '2017-12-08'; 

% Mode: 'Call', 'Put', or 'Both'. 'Both' is best for fixing Skew/Eta.
CalibMode = 'Put'; 

% 2. Locate Target File
filelist = dir(fullfile(data_folder, ['*' target_date '*.csv']));

if isempty(filelist)
    error('No option file found for date %s in %s', target_date, data_folder);
end

% Construct full path for IR_curves
current_file_path = fullfile(filelist(1).folder, filelist(1).name);
fprintf('Processing target file: %s\n', filelist(1).name);

% 3. Load Market Data (IR Curves & Options)
[discounts, fwd, expiries, TradeDate, Spot, Option_data] = IR_curves(current_file_path, spot_file_name);
T_years = years(expiries - TradeDate);

if isempty(T_years)
    error('Data loaded but T_years is empty.');
end

% Initialize marketdata_08 cell (9 columns for hybrid logic)
marketdata_08 = cell(length(Option_data.datesExpiry), 9); 

% 4. Data Preparation & ATM Filtering
fprintf('Applying ATM Delta Filter (0.40 - 0.60)...\n');

for i = 1:length(Option_data.datesExpiry)
    
    Dt = T_years(i);
    B = discounts(i);
    F0 = fwd(i);
    
    if Dt <= 0, continue; end
    
    % --- Extract Raw Data (Force Columns) ---
    RawStrikes = Option_data.strikes(i).value(:);
    RawCallPrices = (Option_data.callAsk(i).prices(:) + Option_data.callBid(i).prices(:))/2;
    RawPutPrices  = (Option_data.putAsk(i).prices(:)  + Option_data.putBid(i).prices(:))/2;
    
    % --- Calculate Deltas (Black-76 Proxy) ---
    sigma_proxy = 0.20;
    d1 = (log(F0 ./ RawStrikes) + (sigma_proxy^2/2)*Dt) / (sigma_proxy*sqrt(Dt));
    
    Delta_Calls = B * normcdf(d1);
    Delta_Puts  = Delta_Calls - B; 
    
    % --- Apply ATM Filter ---
    % Keep Calls where Delta is 0.40 to 0.60
    idx_calls = find(Delta_Calls >= 0.10 & Delta_Calls <= 0.90);
    
    % Keep Puts where |Delta| is 0.40 to 0.60
    idx_puts  = find(abs(Delta_Puts) >= 0.10 & abs(Delta_Puts) <= 0.90);
    
    % --- Combine Data Based on Mode ---
    switch CalibMode
        case 'Call'
            MixedStrikes = RawStrikes(idx_calls);
            MixedPrices  = RawCallPrices(idx_calls);
            TypeFlags    = ones(length(idx_calls), 1); % 1=Call
            
        case 'Put'
            MixedStrikes = RawStrikes(idx_puts);
            MixedPrices  = RawPutPrices(idx_puts);
            TypeFlags    = 2 * ones(length(idx_puts), 1); % 2=Put
            
        case 'Both'
            MixedStrikes = [RawStrikes(idx_calls); RawStrikes(idx_puts)];
            MixedPrices  = [RawCallPrices(idx_calls); RawPutPrices(idx_puts)];
            TypeFlags    = [ones(length(idx_calls), 1); 2 * ones(length(idx_puts), 1)];
    end
    
    % --- Store in Cell Array ---
    marketdata_08{i,1} = Dt; 
    marketdata_08{i,2} = B; 
    marketdata_08{i,3} = F0; 
    marketdata_08{i,4} = MixedPrices; 
    marketdata_08{i,5} = MixedStrikes;
    marketdata_08{i,9} = TypeFlags; % Flag column for Objective Function
end

% 5. Optimization
x0 = [0.15, -1.0, 1.0]; 

fprintf('Starting Calibration (%s Mode)...\n', CalibMode);

% Define Objective Function
objective_function = @(x) Calibration_Objective(x, marketdata_08);

% --- CONSTRAINTS ---
% Format: [sigma, eta, k]
% Lower Bounds: sigma > 1%, eta can be large negative, k > 0
lb = [0.01, -5.0, 0.01]; 

% Upper Bounds: sigma < 100%, ETA MUST BE NEGATIVE (<= -0.01), k < 5
ub = [1.0, -0.01, 5.0]; 

% Setup fmincon options
options = optimoptions('fmincon', 'Display', 'iter', ...
    'Algorithm', 'sqp', ... % SQP is often robust for this
    'DiffMinChange', 1e-4, ...
    'MaxFunctionEvaluations', 3000);

% Run Constrained Optimization
[optimal_params, fval] = fmincon(objective_function, x0, [],[],[],[], lb, ub, [], options);

% 6. Display & Save Results
sigma_opt = optimal_params(1);
eta_opt   = optimal_params(2);
k_opt     = optimal_params(3);

disp('-----------------------------');
disp(['Results for ' datestr(TradeDate)]);
disp(['Spot Price:       ' num2str(Spot)]);
disp(['Calibrated Sigma: ' num2str(sigma_opt)]);
disp(['Calibrated Eta:   ' num2str(eta_opt)]);
disp(['Calibrated k:     ' num2str(k_opt)]);
disp(['Total Error:      ' num2str(fval)]);
disp('-----------------------------');

%%
% --- 7. REGENERATE DENSE GRID FOR PLOTTING ---
% The filtered marketdata_08 has very few points (only ATM). 
% We need to generate a dense grid using the calibrated parameters to plot a surface.

fprintf('Regenerating dense surface data for plotting...\n');

% Settings for Inverse Black-76
tol = 1e-6; 
max_iter = 200;

for i = 1:length(Option_data.datesExpiry)
    Dt = marketdata_08{i,1};
    B  = marketdata_08{i,2};
    F0 = marketdata_08{i,3};
    
    if isempty(Dt) || Dt <= 0
        continue;
    end
    
    % 1. Create a Dense Strike Grid (e.g., 30 points from 70% to 130% of Spot)
    % This ensures linspace inside vol_surface_plot receives a valid vector.
    K_dense = linspace(Spot*0.95, Spot*1.05, 30); 
    
    % 2. Calculate Model Prices (Calls) on this grid
    % We use the calibrated parameters 'optimal_params'
    C_model = NIG_Pricer(optimal_params, F0, K_dense, B, Dt, 1);
    
    % 3. Invert to get Model Implied Volatilities
    iv_model = zeros(size(C_model));
    for k = 1:length(C_model)
        iv_model(k) = Black76Inverse(C_model(k), F0, K_dense(k), B, Dt, tol, max_iter);
    end
    
    % 4. Store in columns 7 (IV) and 8 (Strikes) as expected by vol_surface_plot
    marketdata_08{i, 7} = iv_model; 
    marketdata_08{i, 8} = K_dense;
end


% --- 8. PLOTTING ---
fig = figure('Name', 'NIG Calibrated Volatility Surface', 'Position', [100 100 1200 600]);
main_ax = axes('Parent', fig); 
% Now this will work because Columns 7 and 8 are populated with dense data
vol_surface_plot(marketdata_08, main_ax, TradeDate, Spot);

%%
run_calibration();

%% Vol Surface Animation
vol_surface_evo();

%%
NIG_params = load("./paramsCalibratedParams.mat");
plot_nig_parameter_evo()


%% Question c : Simulation 

% 2. Model & Grid Parameters
% FFT grid parameter (dz)
param = 0.01; 
N_sim=10^7;
% NIG Model Parameters
sigma = sigma_opt ;  % Volatility of the Gaussian component
eta = eta_opt;      % Controls skewness
k = k_opt;        % Controls kurtosis
alpha = 1/2;  % Controls the heaviness of the tails
Dt = 1;       % Time to maturity (in years)
% Market Parameters
F0 = Spot;       % Forward price (S0 if r=q=0)
B = discounts(1)+0.01;        % Discount factor (B = exp(-r*Dt)), implies r=0
%% 3. Model Definition (Characteristic Function)
% Calculate the optimal shift parameter 'a' for the Lewis formula.
% 'p' is the boundary of the analiticity domain.
p = (1/2 + eta) + sqrt((1/2 + eta)^2 + 2 * (1 - alpha) / (sigma^2 * k));
a = p / 2; % Damping parameter
% Define the function L(z) related to the NIG Laplace exponent
if (alpha == 1/2)
    % Standard NIG parameterization
    L = @(z, k, sigma) exp(Dt / k * (1 - sqrt(1 + 2 * k * z * sigma^2)));
else
    % A different parameterization (e.g., related to VG)
    L = @(z, k, sigma) exp(-Dt / k * log(1 + k * z * sigma^2));
end
% Define the risk-neutral Characteristic Function (phi)
phi = @(u) exp(-1i * u * log(L(eta, k, sigma))) .* ...
           L(((u).^2 + 1i * u * (1 + 2 * eta)) / 2, k, sigma);
% --- Sanity Check (Optional) ---
% Check the martingale property: phi(-i) should equal 1
% phi(-1i) 
% ---
%% 4. Setup for Pricing and Error Analysis
% Define the range of strike prices (K) and log-moneyness (kk)
K = linspace(2350, 2700, 40);
kk = log(K ./ F0);
% Initialize arrays to store errors
discr_error = zeros(1, 15 - 4 + 1);
Trunk_error = zeros(1, 15 - 4 + 1);
prices_mat = zeros(15 - 4 + 1, length(K));
%% 5. Error Analysis Loop
% Loop over different values of M (N = 2^M) to analyze convergence
fprintf('Running error analysis for M = 4 to 15...\n');
for m = 4:15
    N = 2^m;
    loop_idx = m - 3;
    
    % 5a. Estimate theoretical errors (for a standard FFT pricing method)
    [discr_error(loop_idx), Trunk_error(loop_idx)] = Error_bound_CDF(...
        sigma, k, eta, alpha, phi, N, 2*pi/param, Dt, 0, 1.5, p);
    % 5b. Simulate random samples from the Characteristic Function (CF)
    sample=SimulateFromCF(phi,m, param, a, N_sim);
    % 5c. Compute MC option prices from the simulated samples
    prices_mat(loop_idx, :) = mean(max(exp(sample)-K,0)) * B * F0;
    MC_SD(loop_idx)=mean(sqrt(var(max(exp(sample)-K,0))/N_sim));
end
fprintf('Error analysis complete.\n');
%% 6. Plot Theoretical Error Bounds (for standard FFT pricing)
figure;
h_fig = gcf;
set(h_fig, 'InvertHardcopy', 'off');
set(h_fig, 'Renderer', 'painters');
hold on;
plot(4:15, log10(discr_error), 'LineWidth', 3, 'Marker', 'o');
plot(4:15, log10(Trunk_error), 'LineWidth', 3, 'Marker', 'x');
hold off;
xlabel('M (N = 2^M)');
ylabel('Log(Error)');
legend('Discretization Error (FFT)', 'Truncation Error (FFT)');
title('Theoretical FFT Pricing Error Bounds');
ylim([-60, 30]);
grid on;
%% 7. Final High-Precision Price Calculation (Benchmark)
% Set M to a high value for the "true" price using FFT-based pricing
M = 15; 
% Compute the final, "accurate" option prices via standard Fourier pricing
prices =PriceCallOption(K,F0,B,phi,M,param,a);
%% 8. Plot Model Implied Volatility (IV Smile)
figure;
% blkimpv computes the Black-Scholes implied volatility
h_fig = gcf;
set(h_fig, 'InvertHardcopy', 'off');
set(h_fig, 'Renderer', 'painters');
plot(K, blkimpv(F0, K, 0, Dt, prices), 'LineWidth', 3);
xlabel('Strike Price (K)');
ylabel('Model Implied Volatility');
title('NIG Model IV Smile (from Benchmark Price)');
grid on;
%% 9. Plot Actual MC Error vs. Theoretical FFT Bounds
% Calculate the actual maximum error of the MC simulation vs. the benchmark price
actual_conv_error = max(abs(prices_mat - prices), [], 2);
figure;
h_fig = gcf;
set(h_fig, 'InvertHardcopy', 'off');
set(h_fig, 'Renderer', 'painters');
hold on;
plot(4:15, log10(discr_error), 'LineWidth', 3, 'LineStyle', '--');
plot(4:15, log10(Trunk_error), 'LineWidth', 3, 'LineStyle', ':');
plot(4:15, log10(actual_conv_error), 'LineWidth', 3, 'Marker', 'o');
plot(4:15,log10(MC_SD), 'LineWidth', 3)
hold off;
xlabel('M (N = 2^M)');
ylabel('Log(Error)');
legend('Discretization (Theory, FFT)', 'Truncation (Theory, FFT)', 'Actual Max MC Price Error', 'MC Standard Deviation');
title('MC Simulation Error vs. Theoretical FFT Bounds');
ylim([-10, 10]);
grid on;


%% PRICE_MIAMI_CERTIFICATE.m

% --- 1. Load Data & Calibration ---
if isfile('./params/NIG_Params_20171208.mat')
    load('./params/NIG_Params_20171208.mat'); 
else
    error('Calibration file not found.');
end


S0 = Spot;
% Dates
FinalDate = datetime('09-Dec-2019');
ValuationDate = datetime('08-Dec-2017');
Notional = 15e6;
% --- 3. Simulation Setup ---
T_final = years(FinalDate - ValuationDate);
r_riskfree = interp1(T_years, -log(discounts)./T_years, T_final, 'spline', 'extrap');
F_2y = interp1(T_years, fwd, T_final, 'spline', 'extrap');
q_div = r_riskfree - log(F_2y / S0) / T_final;

fprintf('Pricing Parameters:\n r = %.4f%%\n q = %.4f%%\n', r_riskfree*100, q_div*100);


% Call the function using the variables returned by your calibration step
[Price_Cert, StdErr] = Price_Miami_Certificate_Function(optimal_params, S0, r_riskfree, q_div, ValuationDate);

% Display
fprintf('\n--------------------------------------\n');
fprintf('MIAMI CROCODILE CERTIFICATE PRICE (DISCRETE BARRIER)\n');
fprintf('--------------------------------------\n');
fprintf('Price:       $ %.2f\n', Price_Cert);
fprintf('%% of Notional: %.2f %%\n', (Price_Cert/Notional)*100);
fprintf('Std Error:   $ %.2f\n', StdErr);
fprintf('--------------------------------------\n');



%% RUN_BACKTEST_STRATEGY.m
% Automates the Hedging Strategy Backtest (Question 6)
% Dates: 08 Dec (Start), 11 Dec, 12 Dec, 13 Dec. 

% --- 1. CONFIGURATION ---
Dates = {'2017-12-08', '2017-12-11', '2017-12-12', '2017-12-13'};
Notional = 15e6; % 15 Million
Pos_Cert = -1;   % Short Position

% Storage
PnL_History = [];
Cash_Account = 0; 

% Portfolio State
Pos_Underlying = 0; 
Pos_Put_Option = 0; 
Current_K = 0; 
Current_T = 0; 

fprintf('Starting Strategy: SELLING VEGA (Short ATM Puts)...\n');

for t = 1:length(Dates)
    CurrentDateStr = Dates{t};
    CurrentDate = datetime(CurrentDateStr);
    
    fprintf('\n--- Date: %s ---\n', CurrentDateStr);
    
    % 1. CALIBRATE & GET DATA
    [optimal_params, S0, r, q, fwd, T_curve, ~] = calibrate_for_date(CurrentDateStr);
    
    % 2. VALUATION (Mark-to-Model)
    [Price_Cert, ~] = Price_Miami_Certificate_Function(optimal_params, S0, r, q, CurrentDate);
    Val_Cert = Pos_Cert * Price_Cert;
    
    % Re-price Existing Hedge Option
    Val_Option = 0;
    Price_Put_Curr = 0;
    
    if Pos_Put_Option ~= 0
        T_remain = years(Current_T - CurrentDate);
        if T_remain > 0.005
            Price_Put_Curr = NIG_Pricer(optimal_params, S0, Current_K, exp(-r*T_remain), T_remain, 2);
            Val_Option = Pos_Put_Option * Price_Put_Curr;
        else
            % Expiry
            Cash_Account = Cash_Account + Pos_Put_Option * max(Current_K - S0, 0);
            Pos_Put_Option = 0;
        end
    end
    
    Val_Stock = Pos_Underlying * S0;
    
    % 3. CALCULATE P&L
    if t > 1
        dt_days = days(CurrentDate - datetime(Dates{t-1}));
        Interest = Cash_Account * (r * dt_days/365);
        Cash_Account = Cash_Account + Interest;
        
        Port_Value = Val_Cert + Val_Stock + Val_Option + Cash_Account;
        Daily_PnL = Port_Value - Prev_Port_Value;
        PnL_History(end+1) = Daily_PnL;
        
        fprintf('   Daily P&L: $ %+.2f\n', Daily_PnL);
    else
        Cash_Account = -Val_Cert; % Premium Received
        Port_Value = 0;
    end
    Prev_Port_Value = Port_Value;
    
    % 4. CALCULATE RISK (GREEKS)
    % A. Certificate Risk (The source of Long Vega)
    [Delta_Cert_Unit, Vega_Cert_Unit] = Calculate_Greeks_NIG(optimal_params, S0, r, q, CurrentDate);
    
    Total_Vega_Risk  = Pos_Cert * Vega_Cert_Unit; % Likely Positive (Long Vega)
    Total_Delta_Risk = Pos_Cert * Delta_Cert_Unit; 
    
    fprintf('   Portfolio Vega Exposure: $ %.2f (Need to Sell)\n', Total_Vega_Risk);
    
    % B. Hedge Instrument (New 3M ATM Put)
    T_new = 0.25; 
    K_new = S0; % ATM for Max Vega
    Target_Exp = CurrentDate + days(90);
    
    [Delta_Put, Vega_Put] = Calculate_Greeks_Vanilla(optimal_params, S0, K_new, r, q, T_new);
    
    % 5. EXECUTE STRATEGY: SELL VEGA
    % We want Net Vega = 0.
    % Pos_Put * Vega_Put + Total_Vega_Risk = 0
    % Pos_Put = -Total_Vega_Risk / Vega_Put
    
    Target_Pos_Put = -Total_Vega_Risk / Vega_Put;
    
    % Since Total_Vega_Risk is likely Positive, Target_Pos_Put will be NEGATIVE.
    % This confirms we are SELLING options.
    
    % 6. DELTA HEDGE (Clean up)
    % Pos_Stock * 1 + Pos_Put * Delta_Put + Total_Delta_Risk = 0
    Target_Pos_Stock = -(Total_Delta_Risk + Target_Pos_Put * Delta_Put);
    
    % 7. TRADING & CASH UPDATE
    
    % Roll Option Position
    if Pos_Put_Option ~= 0 && (K_new ~= Current_K)
        % Buy back old puts (Pay Ask)
        % Approx cost: Mid + Spread
        Cost_Close = -Pos_Put_Option * Price_Put_Curr; 
        Cash_Account = Cash_Account - Cost_Close;
    end
    
    % Sell New Puts (Receive Bid)
    Price_New_Put = NIG_Pricer(optimal_params, S0, K_new, exp(-r*T_new), T_new, 2);
    Premium_Received = -Target_Pos_Put * Price_New_Put; % Neg * Neg = Pos
    Cash_Account = Cash_Account + Premium_Received;
    
    % Adjust Stock
    Trade_Stock = Target_Pos_Stock - Pos_Underlying;
    Cash_Account = Cash_Account - (Trade_Stock * S0);
    
    % Update State
    Pos_Put_Option = Target_Pos_Put;
    Pos_Underlying = Target_Pos_Stock;
    Current_K = K_new;
    Current_T = Target_Exp;
    
    fprintf('   Action: Sold %.0f ATM Puts | Net Stock Position: %.0f\n', ...
        abs(Target_Pos_Put), Pos_Underlying);
end

% --- METRICS ---
Metric = mean(PnL_History.^2);
fprintf('\n========================================\n');
fprintf('STRATEGY PERFORMANCE\n');
fprintf('Metric (Mean Squared P&L): %.2f\n', Metric);
fprintf('Root Mean Square Error:    $ %.2f\n', sqrt(Metric));
fprintf('========================================\n');

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