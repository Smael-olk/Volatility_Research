function plot_term_struct(discounts,R_years,forward_prices,TradeDate,Spot,T_years)
    %PLOT_TERM_STRUCT Summary of this function goes here
    %   Detailed explanation goes here
    figure('Name', 'Term Structure', 'Position', [100 100 1200 600]);
    set(0, 'DefaultFigureColor', [0 0 0]); % Black background
    set(0, 'DefaultAxesColor', [0 0 0]);   % Black axes background
    set(0, 'DefaultAxesXColor', [1 1 1]);  % White X-axis labels/ticks
    set(0, 'DefaultAxesYColor', [1 1 1]);  % White Y-axis labels/ticks
    set(0, 'DefaultTextColor', [1 1 1]);   % white text 
    h_fig = gcf;
    set(h_fig, 'InvertHardcopy', 'off');
    set(h_fig, 'Renderer', 'painters');
    
    orange_color = [251/255 139/255 30/255];
    subplot(1, 2, 1);
    plot(T_years, R_years*100, 'Color', orange_color, 'LineWidth', 2, 'MarkerSize', 6, 'Marker', '^');
    
    title(['Zero rate curve: ', datestr(TradeDate, 'yyyy-mm-dd')]);
    xlabel('Time to Maturity (Years)');
    ylabel('Zero rate r(T)');
    
    % Ensure grid lines are visible on black background (optional, but good practice)
    set(gca, 'GridColor', [0.3 0.3 0.3]); % Dark gray grid for visibility
    grid on;
    axis([0 max(T_years)+0.5 0.5 3]); 
    ylim([0,5])
    
    % 2. Forward Price Curve (Right Subplot)
    subplot(1, 2, 2);
    plot(T_years, forward_prices, 'Color', orange_color, 'LineWidth', 2, 'MarkerSize', 6, 'Marker', 's');
    
    title(['Forward Price Curve: Spot=', num2str(Spot, '%.2f')]);
    xlabel('Time to Maturity (Years)');
    ylabel('Forward Price F(t0,T)');
    
    set(gca, 'GridColor', [0.3 0.3 0.3]); % Dark gray grid
    grid on;
    
    axis([0 max(T_years)+0.5 Spot-10 Spot+100]); 
    ylim([Spot*0.95,Spot*1.05])
    saveas(gcf, "./Plots/Yield2017dec08", 'png');
end

