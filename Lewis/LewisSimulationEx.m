%% Study the
% This script compares two methods for pricing European options under the 
% Normal Inverse Gaussian (NIG) model.
%
% 1. A Monte Carlo (MC) simulation, where random samples are drawn 
%    directly from the characteristic function (CF) using the 'SimulateFromCF' 
%    helper function.
% 2. A high-precision benchmark price calculated using a standard
%    Fourier pricing method (e.g., Lewis/FFT) via the 'PriceCallOption' 
%    helper function.
%
% It performs an analysis by:
% - Calculating the MC price for different simulation parameters (N = 2^M).
% - Plotting the convergence error of the MC price against the benchmark price.
% - Plotting the theoretical error bounds (Discretization, Truncation) 
%   associated with a standard FFT pricing approach for context.
% - Plotting the Monte Carlo standard deviation.
%
% Three helper functions:
%   1. SimulateFromCF(...)
%   2. PriceCallOption(...)
%   3. Error_bound_CDF(...)
% ---
% Author: Michele Azzone
% Date: 12/11/25
% ---
%% 1. Initialization
clear all;
close all;
clc;
%% 2. Model & Grid Parameters
% FFT grid parameter (dz)
param = 0.01; 
N_sim=10^7;
% NIG Model Parameters
sigma = 0.1;  % Volatility of the Gaussian component
eta = 1;      % Controls skewness
k = 1;        % Controls kurtosis
alpha = 1/2;  % Controls the heaviness of the tails
Dt = 1;       % Time to maturity (in years)
% Market Parameters
F0 = 1;       % Forward price (S0 if r=q=0)
B = 1;        % Discount factor (B = exp(-r*Dt)), implies r=0
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
K = linspace(0.8, 1.2, 30);
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
plot(K, blkimpv(F0, K, 0, Dt, prices), 'LineWidth', 3);
xlabel('Strike Price (K)');
ylabel('Model Implied Volatility');
title('NIG Model IV Smile (from Benchmark Price)');
grid on;
%% 9. Plot Actual MC Error vs. Theoretical FFT Bounds
% Calculate the actual maximum error of the MC simulation vs. the benchmark price
actual_conv_error = max(abs(prices_mat - prices), [], 2);
figure;
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