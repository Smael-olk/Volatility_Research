function compare_pricing(DayData, T_index, N_sim, M_grid)
%COMPARE_PRICING_METHODS Compares FFT vs MC prices for a specific maturity.
%   Inputs:
%       DayData: Struct containing calibrated results (Params, Fwds, etc.)
%       T_index: Index of the maturity to test (e.g., 5 for the middle one)
%       N_sim:   Number of MC simulations (e.g., 10^6)
%       M_grid:  Grid resolution for FFT (e.g., 14)

    % --- 1. EXTRACT PARAMETERS ---
    sigma = DayData.Params(1);
    eta   = DayData.Params(2);
    k_val = DayData.Params(3);
    
    Dt = DayData.T_years(T_index);
    F0 = DayData.Fwds(T_index);
    B  = DayData.Discounts(T_index);
    
    fprintf('======================================================\n');
    fprintf('PRICING COMPARISON (T = %.2f years, N_sim = 10^%d)\n', Dt, log10(N_sim));
    fprintf('======================================================\n');

    % --- 2. DEFINE MODEL (Characteristic Function) ---
    alpha = 1/2; 
    p = (0.5 + eta) + sqrt((0.5 + eta)^2 + 2 * (1 - alpha) / (sigma^2 * k_val));
    a = p / 2; % Damping
    param = 0.01; % Integration step

    L = @(z, k, s) exp(Dt / k * (1 - sqrt(1 + 2 * k * z * s^2)));
    phi = @(u) exp(-1i * u * log(L(eta, k_val, sigma))) .* ...
               L(((u).^2 + 1i * u * (1 + 2 * eta)) / 2, k_val, sigma);

    % --- 3. SELECT STRIKES (ITM, ATM, OTM) ---
    % Generate 10 specific strikes from 80% to 120% moneyness
    K_list = linspace(F0 * 0.80, F0 * 1.20, 10)';

    % --- 4. CALCULATE PRICES ---
    
    % A. FFT Benchmark (Semi-Analytical)
    % Using a high M=15 for the "True" reference if M_grid is lower, 
    % or use M_grid to compare "Apple to Apple" implementation.
    % Here we use the user's M_grid for the MC generation, but M=15 for FFT reference.
    P_FFT = PriceCallOption(K_list, F0, B, phi, 15, param, a);
    
    % B. Monte Carlo Simulation
    % Simulate log-returns
    samples = SimulateFromCF(phi, M_grid, param, a, N_sim);
    
    % Price for each strike
    P_MC = zeros(length(K_list), 1);
    S_T = F0 * exp(samples); % Terminal prices
    
    for i = 1:length(K_list)
        Strike = K_list(i);
        payoff = max(S_T - Strike, 0);
        P_MC(i) = mean(payoff) * B;
    end

    % --- 5. DISPLAY RESULTS TABLE ---
    Abs_Error = P_MC - P_FFT;
    Rel_Error_Pct = (Abs_Error ./ P_FFT) * 100;
    
    Tbl = table(K_list, P_FFT, P_MC, Abs_Error, Rel_Error_Pct, ...
        'VariableNames', {'Strike', 'FFT_Price', 'MC_Price', 'Diff_Val', 'Rel_Err_Pct'});
    
    disp(Tbl);
    
    % --- 6. PLOT COMPARISON ---
    figure('Name', 'Price Comparison Check', 'Color', 'black','InvertHardcopy','off');
    BBG_Colors.Orange     = [1.0 0.6 0.0];
    % Subplot 1: Actual Prices
    subplot(2,1,1);
    plot(K_list, P_FFT, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    plot(K_list, P_MC, 'r--x', 'LineWidth', 1.5, 'MarkerSize', 8);
    legend('FFT (Benchmark)', 'Monte Carlo');
    title(sprintf('Option Prices (T=%.2f)', Dt));
    grid on; ylabel('Price');
    
    % Subplot 2: Difference
    subplot(2,1,2);
    bar(K_list, Abs_Error, 'FaceColor', [0.7 0.7 0.7],'Facecolor',BBG_Colors.Orange);
    title('Monte Carlo Noise (MC Price - FFT Price)');
    xlabel('Strike'); ylabel('Price Difference');
    grid on;
end