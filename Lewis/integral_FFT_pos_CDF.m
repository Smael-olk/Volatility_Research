function [I_FFT_all,z]=integral_FFT_pos_CDF(M,param,f)
%Computes the integral of the Fourier Transform of f through 
%the Fast Fourier Transform
%
%INPUTS:
%M        
%param    numerical value; could refer to x1 or dz
%f        function to integrate
%
%OUTPUTS
%I_FFT   value of the integral
%
%DEPENDENCIES:
%interp1


%% Parameters declaration
    dz=param;
    N=2^(M);
    zn=dz*(N-1)*0.5;
    z1=-zn;
    dx=2*pi/(N*dz);
    xn=dx*(N-1)+dx/2;
    x1=dx/2;

%% Compute fft
     
z=z1:dz:zn;

x=x1:dx:xn;

f_j=f(x).*exp(-1i*z1*dx*[0:N-1]);


I_FFT_all=dx*exp(-1i*x1*z).*fft(f_j);



end