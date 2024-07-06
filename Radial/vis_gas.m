function vis=vis_gas(T,M,ro_g)
% Lee-Gonzalez-Eakin correlation (Taken from reservoir engineering handbook
% by Tarek Ahmed)-page 73
%Calculate correlation parameters
x=3.5+(986/T);
y=2.4-0.2*x;
k=((9.4+0.02*M)*T^1.5)/(209+19*M+T);
% Calculate viscosity of the gas
vis=1e-4*k*exp(x*(ro_g/62.4)^y);