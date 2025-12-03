function [sample] = SimulateFromCF(phi, M, param, a, N_sim)
% -------------------------------------------------------------------------
% Purpose:
%   Simulates random samples from a probability distribution given its
%   Characteristic Function (CF).
%
% Method:
%   1. Uses the Lewis inversion formula for positive alpha to numerically compute the
%      Cumulative Distribution Function (CDF) on a grid via FFT.
%   2. Applies Inverse Transform Sampling (X = F_inverse(U)) by interpolating
%      the computed CDF.
%
% -------------------------------------------------------------------------
% INPUTS:
%   phi     : Function handle for the Characteristic Function (e.g., @(u) ...).
%             It's assumed that phi(-1i) = 1 (i.e., E[exp(X)] = 1),
%             which is a common normalization for log-asset prices.
%   M       : Number of points for the Fast Fourier Transform (FFT),
%             typically a power of 2 (e.g., 2^14).
%   param   : Discretization parameter for the FFT grid (controls grid
%             spacing and range).
%   a       : Damping/shift parameter for the inversion integral (the "Lewis shift").
%             Must be chosen such that the damped CF phi(u - 1i*a) is integrable.
%   N_sim   : The desired number of random samples to generate.
%
% OUTPUT:
%   sample  : [N_sim x 1] vector of simulated random samples from the
%             distribution.
% -------------------------------------------------------------------------

%% 1. Define the Integrand for CDF Inversion
% This is the integrand from the Gil-Pelaez/Lewis formula for the
% complementary CDF (P(X > x)).
psi_Lewis = @(u) phi(u - 1i*a) ./ (1i*u + a);

%% 2. Compute the CDF via FFT
% This section relies on a helper function 'integral_FFT_pos_CDF'
% (not provided here) which computes the inversion integral.
%
% We expect this function to return:
%   I_a: The complex-valued integral result on the grid.
%   x_0: The corresponding real-space grid (i.e., the 'x' values).
[I_a, x_0] = integral_FFT_pos_CDF(M, param, @(z) psi_Lewis(z));

% Reconstruct the empirical CDF: F(x) = P(X <= x)
% This is based on the formula: F(x) = 1 - P(X > x),
% where P(X > x) = (exp(-a*x) / pi) * Real[Integral(...)] for x > 0
I = real(I_a);
Grid_emp = 1 - exp(-x_0 * a) / pi .* I;

%% 3. Clean Numerical Artifacts from the CDF
% The numerically computed CDF (Grid_emp) might not be perfectly
% monotonic (i.e., always increasing) due to FFT truncation and
% aliasing errors. We must find the "good" region where it is.

Incr = diff(real(Grid_emp)); % Calculate increments (approx. density)

% Find the first index where the CDF *decreases* (Incr < 0) for x > 0.
% This marks the end of the reliable, positive-x part of the grid.
[~, pos_up] = find((Incr < 0) & (x_0(2:end) > 0));
if isempty(pos_up)
    pos_up = length(Grid_emp) - 1; % No artifact found, use the full grid
end

% Find the last index where the CDF *decreases* for x < 0.
% This marks the start of the reliable, negative-x part of the grid.
[~, pos_down] = find(Incr < 0 & (x_0(2:end) < x_0(pos_up(1))));
if isempty(pos_down)
    pos_down = 0; % No artifact found, start from the beginning
end

% Select the indices corresponding to the monotonic, reliable part of the grid
index = (pos_down(end) + 1):(pos_up(1));

%% 4. Prepare Cleaned Grid for Interpolation
% Extract the "clean" CDF values and their corresponding x-values
Grid_b = Grid_emp(index);
xx = x_0(index);

% Enforce theoretical CDF bounds [0, 1] to correct minor numerical
% over/undershoots (e.g., -1e-16 or 1 + 1e-16).
Grid = max(0, min(Grid_b, 1))';

% For inverse interpolation (mapping F(x) -> x), the 'y' values (Grid)
% must be unique.
[Grid_unique, index_un] = unique(Grid);

% Get the corresponding x-values for the unique CDF points
x_un = xx(index_un);

%% 5. Perform Inverse Transform Sampling
% Generate standard uniform random numbers [0, 1]
u = rand(N_sim, 1);

% Pre-allocate the output sample vector for speed
sample = zeros(N_sim, 1);

% Find indices of uniform numbers that fall *within* the computed CDF range
% [min(Grid_unique), max(Grid_unique)]. We will interpolate these.
ix_used = find((u <= max(Grid_unique)) & (u >= min(Grid_unique)));

% The core inversion: map uniform 'u' (CDF value) back to 'x' (sample value)
% This is the numerical equivalent of sample = F_inverse(u)
sample(ix_used) = interp1(Grid_unique, x_un, u(ix_used));

% Handle extrapolation:
% For 'u' values below our grid's minimum, clip to the smallest x-value.
sample(u < min(Grid_unique)) = min(x_un);

% For 'u' values above our grid's maximum, clip to the largest x-value.
sample(u > max(Grid_unique)) = max(x_un);

end
