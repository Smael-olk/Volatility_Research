function [Price]=PriceCallOption(Strike,F0,B,phi,M,param,a)
%Compute Lewis Formula for a!=0 for European Call option.
%IMPUT:
%Strike   %Strike vector
%F0       %Forward
%B        %Discount
%phi      %Chararcteristic function s.t. phi(-1i)=1
%M        Number of FFT points
%param    Discretization in the real space
%a        Lewis Shift
kk=log(Strike/F0);
psi_Lewis = @(z) 1 ./ ((a + 1i * z) .* (1 + a + 1i * z)) .* ...
                  phi(z - (a + 1) * 1i);
 I = real(integral_FFT_pos(M, param, @(z) psi_Lewis(z), kk));
 Price=(1 * (a < 0) + I .* exp(-a * kk) / (pi)) * B * F0;