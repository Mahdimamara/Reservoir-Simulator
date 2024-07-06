function Bo=oil_fvf(ro_o,gam_oil,Rs,gas_gam)
% See page 97 to find out the uesd correlation to estimate Bo(Advanced Reservoir
% Engineering handbook, Tarek Ahmed)
Bo=(62.4*gam_oil+0.0136*Rs*gas_gam)/(ro_o);