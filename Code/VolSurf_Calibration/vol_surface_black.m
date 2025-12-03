function [s] = vol_surface_black(data, target_ax, now, Spot)
%VOL_SURFACE_PLOT Plots the interpolated Implied Volatility Surface onto specified axes
%   and marks the current Spot Price (T=0) with a projection line to the surface.

    % --- 1. Data Aggregation (unchanged) ---
    all_K = [];
    all_T = [];
    all_IV = [];
    
    num_expiries = size(data, 1);
    
    for i = 1:num_expiries
        T_i = data{i, 1};        
        K_i = data{i, 5};
        IV_i =  data{i, 6};
    
        K_i = K_i(:);  
        IV_i = IV_i(:);
        
        T_vector_i = repmat(T_i, size(K_i, 1), 1); 
        
        all_K = [all_K; K_i];
        all_T = [all_T; T_vector_i];
        all_IV = [all_IV; IV_i];
    end
    
    % --- 2. Grid Creation (unchanged) ---
    K_min = min(all_K);
    K_max = max(all_K);
    T_min = min(all_T);
    T_max = max(all_T);
    
    num_K_points = 300;
    num_T_points = 100;
    
    T_end = max(T_max, 0.01);
    
    K_grid = linspace(K_min, K_max, num_K_points);
    T_grid = linspace(T_min, T_end, num_T_points);
    
    [K_mesh, T_mesh] = meshgrid(K_grid, T_grid);
    
    % --- 3. Interpolation (unchanged) ---
    IV_Surface = griddata(all_K, all_T, all_IV, K_mesh, T_mesh, 'cubic');
      
    % --- 4. Plotting onto Target Axes ---
    cla(target_ax); 

    % Plot the volatility surface (Note: surf doesn't need a legend entry)
    s = surf(target_ax, K_mesh, T_mesh, IV_Surface * 100, 'DisplayName', 'Implied Volatility Surface'); 

    % Get initial Z-limits before plotting the line to use for the marker base
    z_limits = get(target_ax, 'ZLim'); 
    
    % Use 'hold on' for multiple plots
    hold(target_ax, 'on'); 
    
    % --- 5. Add Spot Price Marker on the K-Axis (T=0) ---
    
    X_spot = Spot; 
    Y_spot = 0; 
    Z_spot = 0; % Base of the plot
      
    % Plot the red marker at T=0
    h_spot = plot3(target_ax, X_spot, Y_spot, Z_spot, ...
          'r.', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'Current Spot Price'); 

    % --- 6. Add Projection Line ---
    
    T_proj = linspace(T_min, T_max, 50);
    K_proj = repmat(Spot, size(T_proj));
    IV_proj = griddata(all_K, all_T, all_IV, K_proj, T_proj, 'cubic');
    
    % A. Vertical segment from base (T=0) to first data point (T_min)
    T_vert = [0; T_proj(1)];
    K_vert = [Spot; Spot];
    IV_vert = [Z_spot / 100; IV_proj(1)];

    % Plot the vertical line (Dashed)
    h_vert = plot3(target_ax, K_vert, T_vert, IV_vert * 100, ...
          'r--', 'LineWidth', 1.5, 'DisplayName', 'Spot Projection (Vertical)');
    
    % B. Projection line along the surface (Solid)
    h_proj = plot3(target_ax, K_proj, T_proj, IV_proj * 100, ...
          'w--', 'LineWidth', 3, 'DisplayName', 'ATMF Volatility Term Structure');

    % Use 'hold off' to reset the plotting behavior
    hold(target_ax, 'off'); 

    % --- 7. Customization and Legend ---
    
    set(target_ax.Parent, 'Color', [0 0 0]);
    
    title(target_ax, ['Calibrated IV Surface - Date: ' datestr(now, 'dd-mmm-yyyy') ' | Spot: ' num2str(Spot, '%.2f')], 'Color', 'w');
    xlabel(target_ax, 'Strike Price (K)', 'Color', 'w');
    ylabel(target_ax, 'Time to Maturity (T)', 'Color', 'w');
    zlabel(target_ax, 'Implied Volatility (%)', 'Color', 'w');
    
    % Styling and Appearance
    set(s, 'EdgeColor', 'none'); 
    colormap(target_ax, jet);              
    colorbar(target_ax, 'Color', 'w');
    view(target_ax, 30, 45);               
    grid(target_ax, 'on');

    % Fix Axes limits and include T=0
    Y_limits = get(target_ax, 'YLim');
    ylim(target_ax, [min(0, Y_limits(1)) Y_limits(2)]);
    zlim(target_ax, [min(all_IV) * 100 * 0.9, max(all_IV) * 100 * 1.1]);
    
    % --- Axis Tick Modification (Key Change) ---
    
    % 1. Get current tick locations
    current_ticks = get(target_ax, 'XTick');
    
    % 2. Check if the spot price is already a tick mark (allowing for tolerance)
    if all(abs(current_ticks - Spot) > (K_max - K_min) / 100)
        % If Spot is not close to an existing tick, add it
        new_ticks = sort([current_ticks, Spot]);
    else
        new_ticks = current_ticks;
    end
    
    % 3. Set new ticks
    set(target_ax, 'XTick', new_ticks);
    
    % 4. Format labels to include the Spot value with two decimals
    K_labels = arrayfun(@(k) num2str(k, '%.2f'), new_ticks, 'UniformOutput', false);
    
    % 5. Bold the label corresponding to the Spot Price
    spot_idx = find(abs(new_ticks - Spot) < 1e-6, 1); % Find index of Spot
    if ~isempty(spot_idx)
        K_labels{spot_idx} = ['\bf{' K_labels{spot_idx} '}'];
    end
    
    set(target_ax, 'XTickLabel', K_labels, 'XColor', 'w');

    % --- Legend Implementation (Key Change) ---
    % Only show the spot marker and the projection line
    legend(target_ax, [h_spot, h_proj], 'Location', 'best', 'TextColor', 'w', 'Color', [0.1 0.1 0.1]);

    drawnow;
end