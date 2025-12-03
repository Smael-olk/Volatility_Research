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
    for j=1:length(marketdata{i,4})
        T = marketdata{i,1};
        B= marketdata{i,2};
        F= marketdata{i,3};
        C_target = marketdata{i,4}(j);
        K = marketdata{i,5}(j);
        iv(j)= Black76Inverse(C_target, F, K, B, T, tol, max_iter);
    end
    marketdata{i,6} = iv;
end
fig = figure('Name', ' Black76 Market Implied Volatility Surface', 'Position', [100 100 1200 600]);
main_ax = axes('Parent', fig); % Axes handle for plotting
vol_surface_black(marketdata,main_ax,TradeDate,Spot);
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
    idx_calls = find(Delta_Calls >= 0.40 & Delta_Calls <= 0.60);
    
    % Keep Puts where |Delta| is 0.40 to 0.60
    idx_puts  = find(abs(Delta_Puts) >= 0.40 & abs(Delta_Puts) <= 0.60);
    
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
x0 = [0.2, -0.5, 0.5]; % Initial Guess [sigma, eta, k]
options = optimset('Display', 'iter', 'TolFun', 1e-4, 'MaxFunEvals', 1000);

fprintf('Starting Calibration (%s Mode)...\n', CalibMode);

% Pass marketdata_08 to objective (No regularization needed for single day)
objective_function = @(x) Calibration_Objective(x, marketdata_08);

[optimal_params, fval] = fminsearch(objective_function, x0, options);

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
%%

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
NIG_params = load("CalibratedParams.mat");
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



%%
% --- 2. Product Specifications ---
Notional    = 15e6; % 15 Million USD
S0          = Spot;
Strike      = 1.00 * S0;
Barrier     = 0.90 * S0; % 90%
Trigger     = 1.20 * S0; % 120%
Protection  = 0;
Factor      = 0.90;

% Dates
ValuationDate = datetime('08-Dec-2017');
FinalDate     = datetime('09-Dec-2019');

% Autocall Observation Dates
ObsDates = [datetime('09-Apr-2018'); ...
            datetime('08-Aug-2018'); ...
            datetime('10-Dec-2018'); ...
            datetime('08-Apr-2019'); ...
            datetime('08-Aug-2019'); ...
            datetime('09-Dec-2019')]; % 6th date is Final

% Liquidation Prices (Percentages)
LiqPercs = [1.02; 1.03; 1.05; 1.10; 1.15; 1.20];

% --- 3. Simulation Setup ---
% Estimating r and q
% r: Interpolate zero rate for the product maturity (2 years)
T_final = years(FinalDate - ValuationDate);
% Simple interp from your loaded 'discounts' and 'T_years' vectors
r_riskfree = interp1(T_years, -log(discounts)./T_years, T_final, 'linear', 'extrap');

% q: Dividend Yield. 
% Calibrated Forward F = S0 * exp((r-q)T)  => q = r - ln(F/S0)/T
% We use the 2-year forward from calibration to imply q
F_2y = interp1(T_years, fwd, T_final, 'linear', 'extrap');
q_div = r_riskfree - log(F_2y / S0) / T_final;

fprintf('Pricing Parameters:\n r = %.4f%%\n q = %.4f%%\n', r_riskfree*100, q_div*100);

% Simulation Grid
N_sim = 100000; % 100k simulations
N_days = days(FinalDate - ValuationDate); % Daily steps

fprintf('Simulating %d paths with %d daily steps...\n', N_sim, N_days);

% --- 4. Run Simulation ---
% Note: NIG_Simulate must be in your path
tic;
[S_paths, TimeGrid] = NIG_Simulate(optimal_params, S0, r_riskfree, q_div, T_final, N_days, N_sim);
toc;

%--- 5. Payoff Logic ---

Payoffs = zeros(N_sim, 1);
DiscountFactors = zeros(N_sim, 1);

% Map Observation Dates to Path Indices (Days)
ObsIndices = days(ObsDates - ValuationDate); 
% Ensure indices don't exceed path length (numerical safety)
ObsIndices = min(ObsIndices, size(S_paths, 2)-1); 

for i = 1:N_sim
    Path = S_paths(i, :);
    
    % A. Check Barrier (Continuous / Daily)
    % "Barrier Event" = Reference Value < Trigger AND Below Barrier? 
    % Actually Term sheet says: "Barrier Event has occurred" usually means 
    % touching the barrier level (90%) at ANY time.
    % We check if the path ever dropped below 90%.
    BarrierHit = any(Path < Barrier);
    
    % B. Check Autocall (Dates 1 to 5)
    IsAutocalled = false;
    
    for k = 1:5 
        % Check Spot on Observation Date
        % Index + 1 because Path includes t=0
        S_obs = Path(ObsIndices(k) + 1);
        
        if S_obs >= Trigger
            % Autocall Triggered
            Payoffs(i) = LiqPercs(k) * Notional;
            
            % Discount from THIS date
            T_pay = years(ObsDates(k) - ValuationDate);
            DiscountFactors(i) = exp(-r_riskfree * T_pay);
            
            IsAutocalled = true;
            break; 
        end
    end
    
    if IsAutocalled
        continue; % Done with this path
    end
    
    % C. Final Maturity (Date 6)
    S_final = Path(end);
    T_pay = T_final;
    DiscountFactors(i) = exp(-r_riskfree * T_pay);
    
    % Apply Complex Logic from Term Sheet
    if S_final >= Trigger
        % Scenario: Above Trigger (120%) at Maturity
        if BarrierHit
            % "Reference Value above Trigger Level AND Barrier Event has occurred"
            Payoffs(i) = 1.10 * Notional;
        else
            % "Reference Value above Trigger Level" (Implies No Barrier based on context)
            % (Liquidation price + Additional final amount)
            % Liq Price #6 is 120%. Additional is 23%.
            Payoffs(i) = (1.20 + 0.23) * Notional;
        end
        
    else % S_final < Trigger
        if ~BarrierHit
            % "Below Trigger AND Barrier Event has NOT occurred"
            % Pay Liquidation Price (120%)
            Payoffs(i) = 1.20 * Notional;
        else
            % "Below Trigger AND Barrier Event HAS occurred"
            % Max(Protection, Factor * (S / Strike))
            % Protection = 0%, Factor = 90%
            Val = Factor * (S_final / Strike);
            Payoffs(i) = max(Protection, Val) * Notional;
        end
    end
end

% --- 6. Calculate Price ---
PV = Payoffs .* DiscountFactors;
Price = mean(PV);
StdErr = std(PV) / sqrt(N_sim);

fprintf('\n--------------------------------------\n');
fprintf('MIAMI CROCODILE CERTIFICATE PRICE\n');
fprintf('--------------------------------------\n');
fprintf('Price:       $ %.2f\n', Price);
fprintf('%% of Notional: %.2f %%\n', (Price/Notional)*100);
fprintf('Std Error:   $ %.2f\n', StdErr);
fprintf('--------------------------------------\n');