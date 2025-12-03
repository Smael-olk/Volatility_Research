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

[discounts,forward_prices,expiries,TradeDate,Spot,Option_data]=IR_curves('2017-12-08',spot_file_name); % Discounts and Forward Calibration function
T_years = years(expiries - TradeDate);
R_years = -log(discounts(:)') ./ T_years(:)';
% Plotting
plot_term_struct(discounts,R_years,forward_prices,TradeDate,Spot,T_years)

%% Useful function to visualize Term Structure Evolution for further inspection / refining of hedging strategy.
%term_structure_evo;

%% Market Data 

marketdata = cell(length(Option_data.datesExpiry), 8);
for i=1:length(Option_data.datesExpiry)
    marketdata{i,1} = T_years(i);
    marketdata{i,2} = discounts(i);
    marketdata{i,3} = forward_prices(i);
    marketdata{i,4} = (Option_data.callAsk(i).prices + Option_data.callBid(i).prices ) / 2; % observed market prices
    marketdata{i,5} = Option_data.strikes(i).value;
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
% --- 2. Initial Guess ---
% [sigma, eta, k]
% Start with "safe" parameters that usually imply a valid PDF
x0 = [0.2, -0.5, 0.5];
% --- 3. Optimization Settings ---
options = optimset('Display', 'iter', 'TolFun', 1e-4, 'MaxFunEvals', 1000);
% --- 4. Run Optimization ---
disp('Starting Calibration...');
[optimal_params, fval] = fminsearch(@(x) Calibration_Objective(x, marketdata), x0, options);
% --- 5. Display Results ---
sigma_opt = optimal_params(1);
eta_opt   = optimal_params(2);
k_opt     = optimal_params(3);

disp('-----------------------------');
disp(['Calibrated Sigma: ' num2str(sigma_opt)]);
disp(['Calibrated Eta:   ' num2str(eta_opt)]);
disp(['Calibrated k:     ' num2str(k_opt)]);
disp(['Total Error:     ' num2str(fval)]);


for i=1:length(Option_data.datesExpiry)
    F0= marketdata{i,3};
    B = marketdata{i,2};
    Strikes = marketdata{i,5};
    Dt=marketdata{i,1};
    marketdata{i,6} = NIG_Pricer(optimal_params, F0, Strikes, B, Dt);
end

arbitrage_check(marketdata,optimal_params);

sigma = 0.2; % Initial guess (often 0.2 to 0.4 works well)
tol = 1e-6;
max_iter = 200;
for i=1:length(Option_data.datesExpiry)
    iv= [];
    Dt = marketdata{i,1};
    B= marketdata{i,2};
    F= marketdata{i,3};
    K_market = marketdata{i,5};

    num_points = 20;
    K_min = 2300 ;
    K_max = 2800 ;
    K_dense = linspace(K_min, K_max, num_points);
    marketdata{i,8}=K_dense;
    C=NIG_Pricer(optimal_params, F0, K_dense, B, Dt);
    for j=1:length(C)
        iv(j)= Black76Inverse(C(j), F, K_dense(j), B, Dt, tol, max_iter);
    end
    marketdata{i,7} = iv;
end
fig = figure('Name', 'NIG Implied Volatility Surface', 'Position', [100 100 1200 600]);
main_ax = axes('Parent', fig); % Axes handle for plotting

vol_surface_plot(marketdata,main_ax,TradeDate,Spot);
%%
run_calibration();

%% Vol Surface Animation
vol_surface_evo()

%%
NIG_params = load("CalibratedParams.mat");
plot_nig_parameter_evo()






%% Question c : Simulation 

%% 2. Model & Grid Parameters
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