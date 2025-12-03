function [C] = lewis_pricing(F0 ,Strike,B,Dt)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    % FFT grid parameter (dz)
    param = 0.001; %deltaK
    
    % NIG Model Parameters
    sigma = 0.1;  % Volatility of the Gaussian component
    eta = 1;      % Controls skewness
    k = 1;        % Controls kurtosis
    alpha = 1/2;  % Controls the heaviness of the tails
    % Dt Time to maturity (in years)
    M=15;
    
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
    
    % Define the Lewis formula integrand (psi)
    % This is the dampened characteristic function evaluated on the integration contour
    psi_Lewis2 = @(z) 1 ./ ((a + 1i * z) .* (1 + a + 1i * z)) .* ...
                      phi(z - (a + 1) * 1i);
    
    
    C= PriceCallOption(Strike,F0,B,phi,M,param,a);
end

