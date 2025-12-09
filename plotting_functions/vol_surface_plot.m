function [s] = vol_surface_plot(data, target_ax, now, Spot)
%VOL_SURFACE_PLOT Plots the interpolated Implied Volatility Surface in MONEYNESS.
%   X-Axis: Moneyness (K / Spot)
%   Z-Axis: Implied Volatility (%)
    
    h_fig = gcf;
    set(h_fig, 'InvertHardcopy', 'off');
    
    % --- 1. Data Aggregation (Converted to Moneyness) ---
    all_K = [];
    all_T = [];
    all_IV = [];
    
    num_expiries = size(data, 1);
    
    for i = 1:num_expiries
        T_i = data{i, 1};        
        K_i = data{i, 8};  % Raw Strikes
        IV_i = data{i, 7};
    
        K_i = K_i(:);  
        IV_i = IV_i(:);
        
        % CONVERSION: Normalize Strike by Spot to get Moneyness
        Moneyness_i = K_i ./ Spot; 
        
        T_vector_i = repmat(T_i, size(Moneyness_i, 1), 1); 
        
        all_K = [all_K; Moneyness_i]; % Store Moneyness instead of Strike
        all_T = [all_T; T_vector_i];
        all_IV = [all_IV; IV_i];
    end
    
    % --- 2. Grid Creation ---
    K_min = min(all_K);
    K_max = max(all_K);
    T_min = min(all_T);
    T_max = max(all_T);
    
    num_K_points = 200;
    num_T_points = 100;
    
    T_end = max(T_max, 0.01);
    
    K_grid = linspace(K_min, K_max, num_K_points);
    T_grid = linspace(T_min, T_end, num_T_points);
    
    [K_mesh, T_mesh] = meshgrid(K_grid, T_grid);
    
    % --- 3. Interpolation ---
    IV_Surface = griddata(all_K, all_T, all_IV, K_mesh, T_mesh, 'cubic');
      
    % --- 4. Plotting onto Target Axes ---
    cla(target_ax); 
    s = surf(target_ax, K_mesh, T_mesh, IV_Surface * 100, 'DisplayName', 'Implied Volatility Surface'); 
    
    % Use 'hold on' for multiple plots
    hold(target_ax, 'on'); 
    
    % --- 5. Add Spot Price Marker (Moneyness = 1.0) ---
    
    X_spot = 1.0; % Spot is always 1.0 in Moneyness terms
    Y_spot = 0; 
    Z_spot = 0; 
      
    h_spot = plot3(target_ax, X_spot, Y_spot, Z_spot, ...
          'r.', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'ATMF (100%)'); 
          
    % --- 6. Add Projection Line ---
    
    T_proj = linspace(T_min, T_max, 50);
    K_proj = repmat(1.0, size(T_proj)); % Projection along Moneyness = 1.0
    
    % Interpolate along the 1.0 Moneyness line
    IV_proj = griddata(all_K, all_T, all_IV, K_proj, T_proj, 'cubic');
    
    % A. Vertical segment from base (T=0) to first data point
    T_vert = [0; T_proj(1)];
    K_vert = [1.0; 1.0];
    IV_vert = [Z_spot / 100; IV_proj(1)];
    
    h_vert = plot3(target_ax, K_vert, T_vert, IV_vert * 100, ...
          'r--', 'LineWidth', 1.5, 'DisplayName', 'Spot Projection');
    
    % B. Projection line along the surface
    h_proj = plot3(target_ax, K_proj, T_proj, IV_proj * 100, ...
          'w--', 'LineWidth', 3, 'DisplayName', 'ATMF Term Structure');
          
    hold(target_ax, 'off'); 
    
    % --- 7. Customization and Legend ---
    
    set(target_ax.Parent, 'Color', [0 0 0]);
    
    title(target_ax, ['Calibrated IV Surface - Date: ' datestr(now, 'dd-mmm-yyyy') ' | Spot: ' num2str(Spot, '%.2f')], 'Color', 'w');
    xlabel(target_ax, 'Moneyness (K / Spot)', 'Color', 'w'); % Updated Label
    ylabel(target_ax, 'Time to Maturity (T)', 'Color', 'w');
    zlabel(target_ax, 'Implied Volatility (%)', 'Color', 'w');
    
    % Styling
    set(s, 'EdgeColor', 'none'); 
    colormap(target_ax, jet);              
    colorbar(target_ax, 'Color', 'w');
    view(target_ax, 30, 45);               
    grid(target_ax, 'on');
    
    Y_limits = get(target_ax, 'YLim');
    ylim(target_ax, [min(0, Y_limits(1)) Y_limits(2)]);
    zlim(target_ax, [min(all_IV) * 100 * 0.9, max(all_IV) * 100 * 1.1]);
    
    % --- Axis Tick Modification (Centering on 1.0) ---
    
    current_ticks = get(target_ax, 'XTick');
    
    % Check if 1.0 is in ticks
    if all(abs(current_ticks - 1.0) > 0.01)
        new_ticks = sort([current_ticks, 1.0]);
    else
        new_ticks = current_ticks;
    end
    
    set(target_ax, 'XTick', new_ticks);
    
    % Format labels as percentages (e.g., 0.9 -> 90%, 1.0 -> 100%)
    K_labels = arrayfun(@(k) sprintf('%.0f%%', k*100), new_ticks, 'UniformOutput', false);
    
    % Bold the 100% label
    spot_idx = find(abs(new_ticks - 1.0) < 1e-6, 1);
    if ~isempty(spot_idx)
        K_labels{spot_idx} = ['\bf{' K_labels{spot_idx} '}'];
    end
    
    set(target_ax, 'XTickLabel', K_labels, 'XColor', 'w');
    
    legend(target_ax, [h_spot, h_proj], 'Location', 'best', 'TextColor', 'w', 'Color', [0.1 0.1 0.1]);
    drawnow;
end