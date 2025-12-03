%% 
%% Study the

% This script prices European options using the Normal Inverse Gaussian (NIG)
% model. It uses the Lewis (generalized Fourier) approach, calculating the
% required integral via the Fast Fourier Transform (FFT).
%
% It also performs an error analysis by plotting theoretical error bounds
% (Discretization, Truncation) and the actual convergence error against
% the number of FFT points (N = 2^M).
%
% Assumes the existence of two helper functions:
%   1. integral_FFT_pos(M, param, psi_func, kk, flag)
%   2. Error_bound(sigma, k, eta, alpha, phi_func, N, dx, Dt, r, K, p)

% ---
% Author: Michele Azzone
% Date: 12/11/25
% ---

%% 1. Initialization
clear all;
close all;

%% 2. Model & Grid Parameters

% FFT grid parameter (dz)
param = 0.001; %deltaK

% NIG Model Parameters
sigma = 0.1;  % Volatility of the Gaussian component
eta = 1;      % Controls skewness
k = 1;        % Controls kurtosis
alpha = 1/2;  % Controls the heaviness of the tails
Dt =1;       % Time to maturity (in years)

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

% Define the Lewis formula integrand (psi)
% This is the dampened characteristic function evaluated on the integration contour
psi_Lewis2 = @(z) 1 ./ ((a + 1i * z) .* (1 + a + 1i * z)) .* ...
                  phi(z - (a + 1) * 1i);

%% 4. Setup for Pricing and Error Analysis

% Define the range of strike prices (K) and log-moneyness (kk)
K = linspace(0.6, 1.5, 30);
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
    
    % 5a. Estimate theoretical errors
    [discr_error(loop_idx), Trunk_error(loop_idx)] = Error_bound(...
        sigma, k, eta, alpha, phi, N, 2*pi/param, Dt, 0, 1.5, p);

    % 5b. Compute the FFT integral
    % This implements Eq. (6.1) from the referenced paper.
    I = real(integral_FFT_pos(m, param, @(z) psi_Lewis2(z), kk)); 

    % 5c. Compute option prices using the Lewis formula
    prices_mat(loop_idx, :) = (1 * (a < 0) + I .* exp(-a * kk) / (pi)) * B * F0;
end
fprintf('Error analysis complete.\n');

%% 6. Plot Theoretical Error Bounds

figure;
hold on;
plot(4:15, log10(discr_error), 'LineWidth', 3, 'Marker', 'o');
plot(4:15, log10(Trunk_error), 'LineWidth', 3, 'Marker', 'x');
hold off;
xlabel('M (N = 2^M)');
ylabel('Log(Error)');
legend('Discretization Error', 'Truncation Error');
title('Theoretical FFT Error Bounds');
ylim([-60, 30]);
grid on;

%% 7. Final High-Precision Price Calculation

% Set M to a high value for the "true" price
M = 15; 
param=0.001
% Compute the integral with the final M
I = real(integral_FFT_pos(M, param, @(z) psi_Lewis2(z), kk)); 

% Compute the final, "accurate" option prices
prices = (1 * (a < 0) + I .* exp(-a * kk) / (pi)) * B * F0;

%% 8. Plot Model Implied Volatility (IV Smile)

figure;
% blkimpv computes the Black-Scholes implied volatility
plot(K, blkimpv(F0, K, 0, Dt, prices), 'LineWidth', 3);
xlabel('Strike Price (K)');
ylabel('Model Implied Volatility');
title('NIG Model IV Smile');
grid on;

%% 9. Plot Actual vs. Theoretical Errors

% Calculate the actual maximum error at each M relative to the M=15 price
actual_conv_error = max(abs(prices_mat - prices), [], 2);

figure;
hold on;
plot(4:15, log10(discr_error), 'LineWidth', 3, 'LineStyle', '--');
plot(4:15, log10(Trunk_error), 'LineWidth', 3, 'LineStyle', ':');
plot(4:15, log10(actual_conv_error), 'LineWidth', 3, 'Marker', 'o');
hold off;
xlabel('M (N = 2^M)');
ylabel('Log(Error)');
legend('Discretization (Theory)', 'Truncation (Theory)', 'Actual Max Price Error');
title('Actual Convergence vs. Theoretical Bounds');
ylim([-60, 30]);
grid on;