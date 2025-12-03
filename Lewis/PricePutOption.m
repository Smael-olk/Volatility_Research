function [Price]=PricePutOption(Strike,F0,B,phi,M,param,a)
%Compute Lewis Formula for a!=0 for European Call option.
%IMPUT:
%Strike   %Strike vector
%F0       %Forward
%B        %Discount
%phi      %Chararcteristic function s.t. phi(-1i)=1
%M        Number of FFT points
%param    Discretization in the real space
%a        Lewis Shift
a_put = a + 1;
kk=log(Strike/F0);
psi_Lewis_Put = @(z) 1 ./ ((a + 1i * z) .* (a + 1 + 1i * z)) .* ...
                      phi(z - a_put * 1i);
I = real(integral_FFT_pos(M, param, @(z) psi_Lewis_Put(z), kk));
Price = (1 * (a > -1) - I .* exp(-(a + 1) * kk) / (pi)) * B * F0;
end