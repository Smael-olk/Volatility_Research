function [Price, StdErr] = Price_Miami_Certificate_Function(params, S0, r, q, ValuationDate)
%PRICE_MIAMI_CERTIFICATE_FUNCTION Prices the certificate given model parameters.
%   INPUTS:
%       params        : [sigma, eta, k] (NIG Calibrated Parameters)
%       S0            : Current Spot Price
%       r             : Risk-free rate (constant for simulation horizon)
%       q             : Dividend yield (constant for simulation horizon)
%       ValuationDate : Current date (datetime object)
%
%   OUTPUTS:
%       Price         : Estimated Price of the Certificate
%       StdErr        : Standard Error of the Monte Carlo estimate

    % --- 1. Product Specifications ---
    Notional    = 15e6; 
    Strike      = 1.00 * S0;
    Barrier     = 0.90 * S0; 
    Trigger     = 1.20 * S0; 
    Protection  = 0;
    Factor      = 0.90;
    
    % Dates
    FinalDate = datetime('09-Dec-2019');
    
    % ALL Observation Dates (Autocall + Final Barrier Check)
    ObsDates = [datetime('09-Apr-2018'); ...
                datetime('08-Aug-2018'); ...
                datetime('10-Dec-2018'); ...
                datetime('08-Apr-2019'); ...
                datetime('08-Aug-2019'); ...
                datetime('09-Dec-2019')]; % Date 6 is Final
    
    % Liquidation Prices (Percentages)
    LiqPercs = [1.02; 1.03; 1.05; 1.10; 1.15; 1.20];
    
    % --- 2. Simulation Setup ---
    T_final = years(FinalDate - ValuationDate);
    
    if T_final <= 0
        % Handle expiry case: Payoff is deterministic if we are AT maturity
        % (Assuming we handle this externally, or return 0 for now)
        Price = 0; StdErr = 0; return;
    end

    % Simulation Grid
    N_sim = 100000; 
    N_days = days(FinalDate - ValuationDate); 
    
    % --- 3. Run Simulation ---
    % Generate paths using the NIG model
    [S_paths, ~] = NIG_Simulate(params, S0, r, q, T_final, N_days, N_sim);
    
    % --- 4. Payoff Logic (Discrete Barrier) ---
    Payoffs = zeros(N_sim, 1);
    DiscountFactors = zeros(N_sim, 1);
    
    % Calculate Path Indices for Observation Dates
    % (+1 because S_paths starts at t=0, so Day 1 is index 2)
    ObsIndices = round(days(ObsDates - ValuationDate)) + 1;
    
    % Safety: Ensure indices don't exceed simulation bounds (e.g. leap years)
    ObsIndices = min(ObsIndices, size(S_paths, 2));
    
    for i = 1:N_sim
        Path = S_paths(i, :);
        
        % --- A. CHECK DISCRETE BARRIER ---
        % Only check S_t on the 6 observation dates
        ObservedPrices = Path(ObsIndices);
        BarrierHit = any(ObservedPrices < Barrier);
        
        % --- B. Check Autocall (Dates 1 to 5) ---
        IsAutocalled = false;
        
        for k = 1:5 
            % Check if observation date has passed relative to valuation?
            % Assuming ValuationDate < ObsDates(1) for this logic.
            % Ideally, filter ObsDates for future dates only.
            
            S_obs = Path(ObsIndices(k));
            
            if S_obs >= Trigger
                % Autocall Triggered
                Payoffs(i) = LiqPercs(k) * Notional;
                
                % Discount from THIS date
                T_pay = years(ObsDates(k) - ValuationDate);
                DiscountFactors(i) = exp(-r * T_pay);
                
                IsAutocalled = true;
                break; 
            end
        end
        
        if IsAutocalled
            continue; 
        end
        
        % --- C. Final Maturity (Date 6) ---
        S_final = Path(end);
        T_pay = T_final;
        DiscountFactors(i) = exp(-r * T_pay);
        
        % Apply Maturity Logic
        if S_final >= Trigger
            if BarrierHit
                % Scenario 3: Above Trigger + Barrier Hit -> 110%
                Payoffs(i) = 1.10 * Notional;
            else
                % Scenario 1: Above Trigger + No Barrier -> 120% + 23%
                Payoffs(i) = (1.20 + 0.23) * Notional;
            end
            
        else % S_final < Trigger
            if ~BarrierHit
                % Scenario 2: Below Trigger + No Barrier -> 120%
                Payoffs(i) = 1.20 * Notional;
            else
                % Scenario 4: Below Trigger + Barrier Hit -> Loss
                Val = Factor * (S_final / Strike);
                Payoffs(i) = max(Protection, Val) * Notional;
            end
        end
    end
    
    % --- 5. Calculate Price ---
    PV = Payoffs .* DiscountFactors;
    Price = mean(PV);
    StdErr = std(PV) / sqrt(N_sim);

end