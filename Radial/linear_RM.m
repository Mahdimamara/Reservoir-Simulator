% To generate pressure distribution curves in the reservoir
% Linear flow
clear
clc
fprintf('--------------Linear flow--------------\n')
% Get input data
fprintf('Incompressibile fluid: 1\n');
fprintf('Slightly compressible fluid: 2\n');
fprintf('Compressible fluid: 3\n');
t_f=input('Type of fluid: ');
xf=input('X_final (ft): ');
fprintf('\n');
fprintf('User-input data: 1\n');
fprintf('Correlations: 2\n');
t_i=input('Select type of input data: \n');
% Do calculations
x=linspace(0,xf,1000);
fprintf('Please provide rock and fluid data\n');
fprintf('\n');
switch t_f
    case{1}
        switch t_i
            case{1} 
                k=input('Permeability (md): ');
                vis=input('Viscosity (cp): ');
                A=input('Area (ft^2): ');
                pwf=input('Wellbore pressure (psi): ');
                q=input('Well flow rate (STB/day): ');
                Bo=input('Oil FVF (Rbbl/STB): ');
            case{2}
                k=input('Permeability (md): ');
                A=input('Area (ft^2): ');
                pwf=input('Wellbore pressure (psi): ');
                q=input('Well flow rate (STB/day): ');
                fprintf('Other properties are calculated using appropriate correlations and following inputs\n');
                T=input('Temperature (R): ');
                API=input('Oil API: ');
                vis=vis_oil(T,API);
                fprintf('Oil viscosty was estimated by Beggs and Robinson correlation\n');
                % Calculate Oil FVF using Material balance equation
                % Calculare oil specific gravity
                oil_gam=141.5/(131.5+API);
                % Calculate oil density (pcf)
                ro_o=62.4*oil_gam;
                % Input Rs and gam_gas
                Rs=input('Rs (Scf/STB): ');
                gas_gam=input('Gas specific gravity: ');
                Bo=oil_fvf(ro_o,oil_gam,Rs,gas_gam);
                fprintf('Oil FVF was estimated by Material balance method\n');
        end
        p=((q*Bo*vis*x)/(0.001127*k*A))+pwf;
    case{2}
        switch t_i
            case{1} 
                k=input('Permeability (md): ');
                vis=input('Viscosity (cp): ');
                A=input('Area (ft^2): ');
                co=input('Oil compressibility (1/psi): ');
                Bo=input('Oil FVF (Rbbl/STB): ');
                pwf=input('Wellbore pressure (psi): ');
                q=input('Well flow rate (STB/day): ');
            case{2}
                k=input('Permeability (md): ');
                A=input('Area (ft^2): ');
                pwf=input('Wellbore pressure (psi): ');
                q=input('Well flow rate (STB/day): ');
                fprintf('Other properties are calculated using appropriate correlations and following inputs\n');
                T=input('Temperature (R): ');
                API=input('Oil API: ');
                vis=vis_oil(T,API);
                fprintf('Oil viscosty was estimated by Beggs and Robinson correlation\n');
                % Calculate Oil FVF using Material balance equation
                % Calculare oil specific gravity
                oil_gam=141.5/(131.5+API);
                % Calculate oil density (pcf)
                ro_o=62.4*oil_gam;
                % Input Rs and gam_gas
                Rs=input('Rs (Scf/STB): ');
                gas_gam=input('Gas specific gravity: ');
                Bo=oil_fvf(ro_o,oil_gam,Rs,gas_gam);
                % Calculate co
                % Input Rsb
                Rsb=input('Gas soubility at bubble point (Scf/STB): ');
                co=c_oil(Rsb,gam_gas,API,T,pwf);
                fprintf('Oil compressibility (co) was estimated by Petrosky-Farshad Correlation\n');
        end
        p=pwf+(1-exp((-q*Bo*vis*co*x)/(0.001127*k*A)))/co;
    case{3}
        switch t_i
            case{1} 
                k=input('Permeability (md): ');
                A=input('Area (ft^2): ');
                pwf=input('Wellbore pressure (psi): ');
                q=input('Well flow rate (bbl/day): ');
                vis=input('Gas viscosity (cp): ');
                z=input('Z-factor: ');
                cg=input('Gas isothermal compressibility factor, cg (1/psi): ');
            case{2}
                k=input('Permeability (md): ');
                A=input('Area (ft^2): ');
                pwf=input('Wellbore pressure (psi): ');
                q=input('Well flow rate (bbl/day): ');
                fprintf('Other properties are calculated using appropriate correlations and following inputs\n');
                T=input('Temperature (R): ');
                M=input('Gas MW: ');
                ro_g=input('Gas Density (lbm/ft^3): ');
                vis=vis_gas(T,M,ro_g);
                fprintf('Gas viscosity was estimated using Lee-Gonzalez-Eakin correlation\n');
                Ppc=input('Gas pseudo-critical pressure (psi): ');
                Tpc=input('Gas pseudo-critical temperature (R): ');
                Tpr=T/Tpc;
                % Try and error procedure to calcualte fluid properties
                p2g=pwf+500;
                % Define Standard conditions
                psc=14.7;Tsc=520;
                eps=1;
                while eps>5e-2
                    p_avg=sqrt((pwf^2+p2g^2)/2);
                    Ppr=p_avg/Ppc;
                    z=compr(Tpr,Ppr);
                    cg=c_gas(Tpr,Ppr,p_avg);
                    % Calculate p2
                    p2=sqrt(pwf^2+((q*psc*z*vis*T*xf)/(Tsc*0.003164*A*k)));
                    eps=abs(p2-p2g);
                    p2g=p2;
                end
                fprintf('Gas Z-factor was calculated using Drunchuk-Purvis-Robinson correlation\n');
                fprintf('cg was calculated using reduced gas compressibility concept and with aid of Drunchuk-Purvis-Robinson correlation\n');
        end
        q=q*5.615;        
        p=sqrt(pwf^2+((q*psc*z*vis*T*x)/(Tsc*0.003164*A*k)));
end
% Plot p vs. x
figure
plot(x,p,'--r');
xlabel('x (ft)');
ylabel('pressure (psi)');
qs=input('Do you want to calculate flow rate using pwf and p_res? (Y=1 ;N=2) ');
if qs==1
    pres=input('Reservoir pressure (psi): ');
    pwf=input('Wellbore presure (psi): ');
    switch t_f
        case{1}
            q_c=((pres-pwf)*0.001127*k*A)/(vis*Bo*xf);
            fprintf('q (STBD)= %f\n',q_c);
        case{2}
            q_c=(-0.001127*k*A*ln(1+co*(pwf-pres)))/(vis*Bo*co*xf);
            fprintf('q (STBD)= %f\n',q_c);
        case{3}
            q_c=(0.003164*Tsc*k*A(pres^2-pwf^2))/(psc*T*xf*vis*z);
            fprintf('q (SCF/day)= %f\n',q_c);
    end
else
    ...
end