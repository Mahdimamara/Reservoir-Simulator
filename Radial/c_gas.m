function cg=c_gas(Tpr,Ppr,p)
% See pages 63-65 to find out the uesd correlation to estimate cg(Advanced Reservoir
% Engineering handbook, Tarek Ahmed)
% Define Coefficients
A1=0.31506237;A2=-1.0467099;A3=-0.57832720;A4=0.53530771;A5=-0.61232032;
A6=-0.10488813;A7=0.68157001;A8=0.68446549;
T1=A1+A2/Tpr+A3/(Tpr^3);
T2=A4+A5/Tpr;
T3=(A5*A6)/(Tpr);
T4=A7/(Tpr^3);
T5=(0.27*Ppr)/Tpr;
% Find Ro_r usnig Drunchuk-Purvis-Robinson Correlation
Ro_ri=(0.27*Ppr)/Tpr;
f=@(Ro_r) 1+T1*Ro_r+T2*Ro_r^2+T3*Ro_r^5+(T4*Ro_r^(2)*(1+A8*Ro_r^(2))*exp(-A8*Ro_r^(2)))-T5/Ro_r;
Ro_rans=fzero(f,Ro_ri);
z=(0.27*Ppr)/(Ro_rans*Tpr);
% Calculate gas isothermal compressibility factor
dz_dro=T1+2*T2*Ro_rans+5*T3*Ro_rans^(4)+2*T4*Ro_rans*(1+A8*Ro_rans^(2)-A8^(2)*Ro_rans^(4))*exp(-A8*Ro_rans^(2));
cpr=1/Ppr-((0.27*dz_dro)/(z^(2)*Tpr*(1+(Ro_rans*dz_dro)/z)));
Ppc=p/Ppr;
cg=cpr/Ppc;