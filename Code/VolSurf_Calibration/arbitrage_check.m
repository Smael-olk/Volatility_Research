function arbitrage_check(data,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    figure; 
    hold on;
    title('Risk Neutral Density (RND) Check Observed Market Prices');
    xlabel('Strike (K)');
    ylabel('Density (2nd Derivative)');
    set(gcf, 'InvertHardcopy', 'off');
    set(gcf, 'Renderer', 'painters');
    
    % Loop through each maturity in your data
    num_expiries = size(data, 1);
    
    for i = 1:num_expiries
        % --- 1. Extract Raw Data for this Maturity ---
        Dt = data{i, 1};       % Time to Maturity
        F0 = data{i, 3};       % Forward Price (Check your column index)
        B = data{i, 2};       % Discount Factor (Check your column index)
        K_market = data{i, 5};   % Raw Strikes
        
        % --- 2. Create a Dense, Uniform Strike Grid ---
        % We need this for accurate numerical differentiation
        num_points = 200;
        K_min = min(K_market) * 0.8;
        K_max = max(K_market) * 1.2;

        
        % Create a vector of 200 evenly spaced strikes
        K_dense = linspace(K_min, K_max, num_points);
        dK = K_dense(2) - K_dense(1); % The step size
        
        % --- 3. Interpolate IV onto the Dense Grid ---
        % Use 'pchip' (Piecewise Cubic Hermite) to avoid oscillations that 'spline' 
        % might create. Oscillations cause fake arbitrage!
        C_dense = NIG_Pricer(params, F0, K_dense, B, Dt);
        
        
        % --- 5. Calculate Second Derivative (Finite Difference) ---
        % Formula: (C(k+1) - 2*C(k) + C(k-1)) / dK^2
        % This is the numerical approximation of d^2C/dK^2
        
        RND = diff(C_dense, 2) / (dK^2);
        
        % Adjust K axis for plotting (diff reduces array size by 2)
        K_plot = K_dense(2:end-1);
        
        % --- 6. Plot the RND Slice ---
        plot(K_plot, RND, 'LineWidth', 2, 'DisplayName', ['T = ' num2str(Dt)]);
        
        % --- 7. Check for Arbitrage ---
        if any(RND < 0)
            disp(['WARNING: Arbitrage detected at T = ' num2str(Dt) ...
                  '. Density dips below zero.']);
        end
    end
    
    yline(0, 'r--', 'Zero Bound'); % Add a red line at 0
    legend show;
    grid on;
    hold off;
end

