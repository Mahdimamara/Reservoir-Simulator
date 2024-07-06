function vis=vis_oil(T,API)
% Beggs and Robinson correlation (Taken from reservoir engineering handbook
% by Tarek Ahmed)-page 110
% Calculate Correlation parameters
z=3.0324-0.02023*API;
y=10^(z);
x=y*(T-460)^(-1.163);
% Calculate viscosity
vis=10^(x)-1;