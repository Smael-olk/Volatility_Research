function arbitrage_check(data, params)
%ARBITRAGE_CHECK Checks for RND positivity and calculates Skewness.
    
    figure('Name', 'Risk Neutral Density & Skew Check','InvertHardcopy','off'); 
    hold on;
    title('Risk Neutral Density (PDF) & Skewness');
    xlabel('Strike (K)');
    ylabel('Probability Density');
    
    num_expiries = size(data, 1);
    colors = parula(num_expiries); % Generate distinct colors
    
    fprintf('\n--- SKEWNESS ANALYSIS ---\n');
    
    for i = 1:num_expiries
        % --- 1. Extract Data ---
        Dt = data{i, 1};      
        B  = data{i, 2};      
        F0 = data{i, 3};      
        K_market = data{i, 5}; % Raw Strikes
        
        if isempty(Dt) || Dt <= 0, continue; end
        
        % --- 2. Dense Grid for Derivatives ---
        % Create a wider grid to capture the full tails for accurate moment calc
        K_min = min(K_market) * 0.7;
        K_max = max(K_market) * 1.3;
        num_points = 500; % Higher resolution
        K_dense = linspace(K_min, K_max, num_points);
        dK = K_dense(2) - K_dense(1);
        
        % --- 3. Calculate Prices ---
        % Ensure NIG_Pricer handles the flag 1 (Call)
        C_dense = NIG_Pricer(params, F0, K_dense, B, Dt, 1);
        
        % --- 4. Calculate RND (PDF) ---
        % RND = d^2C / dK^2 * B^-1 (To remove discount scaling if desired, but shape is same)
        % We leave B in to match standard pricing units, or remove it to get true PDF
        PDF = diff(C_dense, 2) / (dK^2) / B; 
        
        % Align K axis (diff reduces size by 2)
        K_pdf = K_dense(2:end-1);
        
        % --- 5. CALCULATE SKEWNESS (3rd Moment) ---
        % Normalize PDF so it sums to 1 (Approximate Integration)
        Prob = PDF * dK;
        Total_Prob = sum(Prob);
        Prob = Prob / Total_Prob; % Force sum to 1 for moment calc
        
        % Mean (Should be close to F0)
        Mean_K = sum(K_pdf .* Prob);
        
        % Variance (2nd Moment)
        Var_K = sum(((K_pdf - Mean_K).^2) .* Prob);
        Std_K = sqrt(Var_K);
        
        % Skewness (3rd Moment)
        % Skew = E[(x - mu)^3] / sigma^3
        Skew_Val = sum(((K_pdf - Mean_K).^3) .* Prob) / (Std_K^3);
        
        % Interpret
        if Skew_Val < -0.1
            SkewType = 'Left (Crash)';
        elseif Skew_Val > 0.1
            SkewType = 'Right (Rally)';
        else
            SkewType = 'Symmetric';
        end
        
        fprintf('T=%.2f yr | Fwd=%.2f | Skew Coeff: %.4f -> %s\n', Dt, F0, Skew_Val, SkewType);
        
        % --- 6. Plotting ---
        plot(K_pdf, PDF, 'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', ['T=' num2str(Dt, '%.2f') ' (' SkewType ')']);
         
        % Add marker for Forward Price on the curve
        % Find PDF value at Forward
        [~, idx_fwd] = min(abs(K_pdf - F0));
        plot(F0, PDF(idx_fwd), 'o', 'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
    end
    
    
    legend('Location', 'best');
    grid on;
    hold off;
end