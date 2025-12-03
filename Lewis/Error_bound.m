function [discr_error,trunc_error] = Error_bound(sigma,k,eta,alpha,f,N,L,t,s,Max_strike,p_a)
%Compute discretization e truncation for alpha>0 for optimal a
%%INPUT
%ATS_params       ATS parameters
%l                ATS laplace trasform 
%N                Discretization parameter
%L                Truncation Param
%D                Regularization Param
%N_sim            Number of simulation
%t                End of increment
%s                Beginning of increment

%%OUTPUT
%errors           Regularization, Truncation and Discretization errors
p=p_a-1;
delta=max(1e-12,L/(N));
N=max(2,N);


%Truncation
    B=exp((1-alpha)/alpha*(t/k-s/k));
    b=(1-alpha)/(alpha*k^(1-alpha))*(sigma^2/(2*(1-alpha)))^alpha...
        *(t^(1)-s^(1));
    omega=2*alpha;
    trunc_error=exp(p/2*Max_strike)*quadgk(@(u)B*exp(-u.^omega*b)./u.^2,L,inf);
%Discretization
discr_error=exp(-2*pi*p/2/delta)+exp(p*Max_strike)*exp(-2*pi*p/2/delta)*f(-1i*(p+0.999));



