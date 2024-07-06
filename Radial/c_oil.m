function co=c_oil(Rsb,gam_gas,API,T,p)
% See page 100 to find out the uesd correlation to estimate co(Advanced Reservoir
% Engineering handbook, Tarek Ahmed)
% The Petrosky-Farshad Correlation
co=1.705e-7*Rsb^(0.69357)*gam_gas^(0.1885)*API^(0.3272)*(T-460)^(0.6729)*p^-0.5906;