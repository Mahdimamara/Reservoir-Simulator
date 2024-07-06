clc
clear all
close all
format bank
%-------------------------------Input-------------------------------------%
%=========================================================================%
Poi=input('Initial Oil Pressure(Psi)[3000]: ');
if isempty(Poi)
    Poi=3000
end
Swi=input('Initial Water Saturation[0.2]: ');
if isempty(Swi)
    Swi=0.2
end
Sor=input('Residual Oil Saturation[0.2]: ');
if isempty(Sor)
    Sor=0.2
end
Phi=input('Porosity[0.25]: ');
if isempty(Phi)
    Phi=0.25
end
K=input('Permeability(Darcy)[0.1]: ');
if isempty(K)
    K=0.1
end
Re=input('Outer radius of Reservoir(ft)(Square)[500]: ');
if isempty(Re)
    Re=500
end
Rw=input('Inner radius of Reservoir(ft)(Square)[2]: ');
if isempty(Rw)
    Rw=2
end
deltaz=input('Reservoir Thickness(ft)[100]: ');
if isempty(deltaz)
    deltaz=100
end
gr=input('Number of Grids in r-direction[5]: ');
if isempty(gr)
    gr=5
end
gte=input('Number of Grids in theta-direction: ');
if isempty(gte)
    gte=6
end
format long
Cr=input('Rock Compressibility[2*10^-6]: ');
if isempty(Cr)
    Cr=2*10^-6
end
format bank
TTD=input('Number Of Changes Of Rates[1]: ');
if isempty(TTD)
    TTD=1
end
disp('HINT: Time Steps Are 1 Hour.');
NT=input('Number of 1st Steps(Hour)[24]: ');
if isempty(NT)
    NT=24
end
Qscoo=input('Oil Rate @ 1st Steps(STBD)[100]: ');
if isempty(Qscoo)
    Qscoo=100
end
for i=2:TTD
NT(i)=input('Number Of Next Steps(Hour): ')
Qscoo(i)=input('Oil Rate @ These Steps(STBD): ')
end
tm=ceil(26.*Qscoo(1)./100);
Qsco=Qscoo./gte;
API=input('Oil Gravity(API) [32]: ')
if isempty(API)
    API=32
end
Gamao=141.5/(API+131.5);
Gamag=input('Specific Gravity of the Solution Gas(gama g) [0.8]: ')
if isempty(Gamag)
   Gamag=0.8
end
TR=input('Temperature of Reservoir(R) [200]: ')
if isempty(TR)
   TR=200
end
%-------------------------Radial grid-------------------------------%
%=========================================================================%
rcoef=(Re/Rw)^(1/gr);
rd=zeros(gr+1,gte);
for j=1:gte
rd(1,j)=Rw;
end
for i=1:gr
    for j=1:gte
    rd(i+1,j)=rcoef*rd(i,j);
    end
end

ud=log(rd);
del_ud=ud(2,1)-ud(1,1);

for j=1:gte
cud(1,j)=0.345;
end
for i=2:gr+1
    for j=1:gte
        cud(i,j)=(ud(i-1,j)+ud(i,j))/2;
    end 
end
crd=exp(cud);
%--------------------------Cell Dimensions--------------------------------%
%=========================================================================%
del_th=360/gte;
%-----------------------------Primary Matrix------------------------------%
%=========================================================================%
Po=Poi.*ones(gr+1,gte);
Sw=(Swi+10^(-10))*ones(gr+1,gte);
Sg=zeros(gr+1,gte);
So=(1-(Swi+10^(-10)))*ones(gr+1,gte);
Bo=ones(gr+1,gte);
Bg=ones(gr+1,gte);
Rs=ones(gr+1,gte);
TroPlus=ones(gr+1,gte);
TthoPlus=ones(gr+1,gte);
TroMinus=ones(gr+1,gte);
TthoMinus=ones(gr+1,gte);
TrgPlus=ones(gr+1,gte);
TthgPlus=ones(gr+1,gte);
TrgMinus=ones(gr+1,gte);
TthgMinus=ones(gr+1,gte);
Cpoo=ones(gr+1,gte);
Cpog=ones(gr+1,gte);
aalpha=ones(gr+1,gte);
l=0;
%----------------------------Calculations---------------------------------%
%=========================================================================%
for SSS=1:TTD
    T=(1/240)*ones(1,NT(SSS));
for iii=1:NT(SSS)
    t(iii)=T(iii);
            Poo=Po;
            %Gas solubility standing correlation
            x=0.0125*API-0.0009*(TR);
            Rs=Gamag.*((Po./18.2+1.4).*(10.^x)).^1.2048;

            RstoPo=(753.*10.^(API./80 -(9.*TR)./10000 ).*Gamag.*...
                (10.^(API./80 -(9.*TR)./10000).*((5.*Po)./91 +7./5)).^(128./625))./11375;

            %Bubble-Point Pressure standing correlation
            a=0.00091.*(TR)-0.0125.*API;
            Pb=18.2.*(((Rs./Gamag).^0.83).*(10.^a)-1.4);

            %Oil Formation Volume Factor standing correlation
            Bo=0.9759+0.00012.*(Rs.*(Gamag./Gamao).^0.5+1.25.*(TR)).^1.2;
            %1/Bo to Po
            
            BotoPo=0.00017;
            Sww=Sw;
            Sgg=Sg;
            %VISCOSITY OF THE DEAD OIL  Beggs-Robinson Correlation
            zz=3.0324-0.02023.*API;
            xx=(10.^zz).*(TR).^-1.163;
            Muod=10.^(xx)-1;
            %VISCOSITY OF THE Saturated OIL  Beggs-Robinson Correlation
            aa=10.715.*(Rs+100).^-0.515;
            bb=5.44.*(Rs+150).^-0.338;
            Muo=aa.*(Muod).^bb;

            for i=1:gr+1
                for j=1:gte
                    if So(i,j)>=0.23 && So(i,j)<=0.8
                        Pc(i,j)=(0.001.*Sg(i,j)).^2;
                        PctoSg(i,j)=0.002.*Sg(i,j);
                        Kro(i,j)=(1-Sg(i,j)/(1-Swi))^4;
                        Krg(i,j)=(Sg(i,j)/(1-Swi))*(2-(Sg(i,j)/(1-Swi)));
                    elseif So(i,j)<=0.23
                        Pc(i,j)=0;
                        Kro(i,j)=0;
                        Krg(i,j)=(Sg(i,j)/(1-Swi))*(2-(Sg(i,j)/(1-Swi)));
                    elseif So(i,j)>=0.8
                        Pc(i,j)=0;
                        Krg(i,j)=0;
                        Kro(i,j)=(1-Sg(i,j)/(1-Swi))^4;
                    end
                        Kro(i,j)=1-Sg(i,j)./20;
                        Krg(i,j)=((Sg(i,j)./(1-Swi)).*(2+(Sg(i,j)./(1-Swi)))./10)^(0.6);
                end
            end 

            
            Pg=Po+Pc;
            Zc=0.85;
            Bg=0.005035.*(TR+460).*Zc./Pg;
            BgtoPg=1./(0.005035.*(TR+460).*Zc);
            %THE VISCOSITY OF NATURAL GASES  Lee-Gonzalez-Eakin Method
            %Assume Methan gas
            Ma=16;
            deng=Pg.*Ma./(10.73.*(TR+460).*Zc);
            xmi=3.5+986./(TR+460)+0.01.*Ma;
            ymi=2.4-0.2.*xmi;
            kmi=((9.4+0.02.*Ma).*(TR+460).^1.5)/(209+19.*Ma+(TR+460));
            Mug=(10.^-4).*kmi.*exp(xmi.*(deng./62.4).^ymi);

            Landao=Kro./((Muo.*Bo));
            Landag=Krg./((Mug.*Bg));
            Vb=zeros(gr+1,gte);
            for j=1:gte
               Vb(1,j)=pi*(Rw^2)*deltaz/gte; 
            end
            for i=2:gr+1
                for j=1:gte
                Vb(i,j)=pi*((rd(i))^2-(rd(i-1))^2)*deltaz/gte;
                end
            end
%--------------------------Transmissibility Terms-------------------------%
%=========================================================================%
   %T theta oPlus--------------------
            for i=1:gr+1
                for j=1:gte-1
                    if Po(i,j+1)>=Po(i,j)
                        TteoPlus(i,j)=2*1.127*Vb(i,j)*K*Landao(i,j+1)/((del_th)^2);
                        RstPlus(i,j)=Rs(i,j+1);
                    else
                        TteoPlus(i,j)=1.127*Vb(i,j)*K*Landao(i,j)/((del_th)^2);
                        RstPlus(i,j)=Rs(i,j);
                    end
                end
            end
             for i=1:gr+1
                    if Po(i,gte)>=Po(i,1)
                        TteoPlus(i,gte)=2*1.127*Vb(i,gte)*K*Landao(i,gte)/((del_th)^2);
                        RstPlus(i,gte)=Rs(i,gte);
                    else
                        TteoPlus(i,gte)=1.127*Vb(i,gte)*K*Landao(i,1)/((del_th)^2);
                        RstPlus(i,gte)=Rs(i,1);
                    end
             end
      
   %T theta oMinus--------------------
            for i=1:gr+1
                for j=2:gte
                    if Po(i,j-1)>=Po(i,j)
                       TteoMinus(i,j)=2*1.127*Vb(i,j)*K*Landao(i,j-1)/((del_th)^2);
                       RstMinus(i,j)=Rs(i,j-1);
                    else
                        TteoMinus(i,j)=2*1.127*Vb(i,j)*K*Landao(i,j)/((del_th)^2); 
                        RstMinus(i,j)=Rs(i,j);
                    end
                end
            end
           for i=1:gr+1
                    if Po(i,gte)>=Po(i,1)
                       TteoMinus(i,1)=2*1.127*Vb(i,1)*K*Landao(i,gte)/((del_th)^2);
                       RstMinus(i,1)=Rs(i,gte);
                    else
                        TteoMinus(i,1)=2*1.127*Vb(i,1)*K*Landao(i,1)/((del_th)^2); 
                        RstMinus(i,1)=Rs(i,1);
                    end
           end
           
   %T r oPlus--------------------
            for i=1:gr
                for j=1:gte
                    if Po(i+1,j)>=Po(i,j)
                       TroPlus(i,j)=1.127*Vb(i,j)*K*Landao(i+1,j)/((del_ud)^2);
                       RsrPlus(i,j)=Rs(i+1,j);
                    else
                        TroPlus(i,j)=1.127*Vb(i,j)*K*Landao(i,j)/((del_ud)^2);
                        RsrPlus(i,j)=Rs(i,j);
                    end
                end
            end
            TroPlus(gr+1,:)=0;
            RsrPlus(gr+1,:)=0;
            for i=2:gr+1
                for j=1:gte
                    if Po(i-1,j)>=Po(i,j)
                       TroMinus(i,j)=1.127*Vb(i,j)*K*Landao(i-1,j)/((del_ud)^2);
                       RsrMinus(i,j)=Rs(i-1,j);
                    else
                        TroMinus(i,j)=1.127*Vb(i,j)*K*Landao(i,j)/((del_ud)^2); 
                        RsrMinus(i,j)=Rs(i,j);
                    end
                end
            end
            TroMinus(1,:)=0;
            RsrMinus(1,:)=Rs(1,:);
   %T theta gPlus--------------------
            for i=1:gr+1
                for j=1:gte-1
                    if Pg(i,j+1)>=Pg(i,j)
                       TtegPlus(i,j)=2*1.127*Vb(i,j)*K*Landag(i,j+1)/((del_th)^2);
                    else
                        TtegPlus(i,j)=2*1.127*Vb(i,j)*K*Landag(i,j)/((del_th)^2); 
                    end
                end
            end
            for i=1:gr+1
                    if Po(i,gte)>=Po(i,1)
                        TtegPlus(i,gte)=2*1.127*Vb(i,gte)*K*Landag(i,gte)/((del_th)^2);
                    else
                        TtegPlus(i,gte)=1.127*Vb(i,gte)*K*Landag(i,1)/((del_th)^2);
                    end
             end
            
   %T theta gMinus--------------------
            for i=1:gr+1
                for j=2:gte
                    if Pg(i,j-1)>=Pg(i,j)
                       TtegMinus(i,j)=2*1.127*Vb(i,j)*K*Landag(i,j-1)/((del_th)^2);
                    else
                        TtegMinus(i,j)=2*1.127*Vb(i,j)*K*Landag(i,j)/((del_th)^2);
                    end
                end
            end
           for i=1:gr+1
                    if Po(i,gte)>=Po(i,1)
                       TtegMinus(i,1)=2*1.127*Vb(i,1)*K*Landag(i,gte)/((del_th)^2);
                    else
                        TtegMinus(i,1)=2*1.127*Vb(i,1)*K*Landag(i,1)/((del_th)^2); 
                    end
           end
            
   %T r gPlus--------------------
            for i=1:gr
                for j=1:gte
                    if Pg(i+1,j)>=Pg(i,j)
                       TrgPlus(i,j)=1.127*Vb(i,j)*K*Landag(i+1,j)/((del_ud)^2);
                    else
                        TrgPlus(i,j)=1.127*Vb(i,j)*K*Landag(i,j)/((del_ud)^2); 
                    end
                end
            end
            TrgPlus(gr+1,:)=0;
   %T r gMinus--------------------
            for i=2:gr+1
                for j=1:gte
                    if Pg(i-1,j)>=Pg(i,j)
                       TrgMinus(i,j)=1.127*Vb(i,j)*K*Landag(i-1,j)/((del_ud)^2);
                    else
                        TrgMinus(i,j)=1.127*Vb(i,j)*K*Landag(i,j)/((del_ud)^2); 
                    end
                end
            end
            TrgMinus(1,:)=0;
%-------------------------------------------------------------------------%

           Cpoo=(Phi.*(1-Sg-Sw).*Vb.*((Cr./Bo)+BotoPo))./(5.615*t(iii));
           Cpog=Phi.*Vb.*(Sg.*((Cr./Bg)+BgtoPg)+Rs.*(1-Sg).*((Cr./Bg)+BgtoPg)...
               +((1-Sg)./Bo).*RstoPo)./(5.615*t(iii));
           Csgg=Phi.*Vb.*(Sg.*((Cr./Bg)+BgtoPg).*(PctoSg)-(Rs./Bo)+1./Bg)./(5.615*t(iii));
           Csgo=-Vb.*Phi./(5.615.*Bo.*t(iii));
           aalpha=-Csgo./Csgg;
           
%------------------------Matirix_Pressure---------------------------------%
%=========================================================================%
           A=TroMinus+aalpha.*(TrgMinus+RsrMinus.*TroMinus);
           E=TteoMinus+aalpha.*(TtegMinus+RstMinus.*TteoMinus);
           C=TroPlus+aalpha.*(TrgPlus+RsrPlus.*TroPlus);
           D=TteoPlus+aalpha.*(TtegPlus+RstPlus.*TteoPlus);
           B=-(TroMinus+TteoMinus+TroPlus+TteoPlus+Cpoo)-aalpha.*((TrgMinus+RsrMinus.*TroMinus)+...
              (TtegMinus+RstMinus.*TteoMinus)+(TrgPlus+RsrPlus.*TroPlus)+...
              (TtegPlus+RstPlus.*TteoPlus)+Cpog);
           
           for j=1:gte
               Pc(1,j)=0;
           end
                      
           for i=2:gr
               for j=2:gte-1
               F(i,j)=-(aalpha(i,j)*(TrgPlus(i,j)*(Pc(i+1,j)-Pc(i,j))+...
               TrgMinus(i,j)*(Pc(i-1,j)-Pc(i,j))+TtegPlus(i,j)*(Pc(i,j+1)-Pc(i,j))...
               +TtegMinus(i,j)*(Pc(i,j-1)-Pc(i,j))))...
               -Po(i,j)*(Cpoo(i,j)+aalpha(i,j)*Cpog(i,j));
               end
           end
           
               for j=2:gte-1
               F(1,j)=-(aalpha(1,j)*(TrgPlus(1,j)*(Pc(2,j)-Pc(1,j))+...
               TtegPlus(1,j)*(Pc(1,j+1)-Pc(1,j))...
               +TtegMinus(1,j)*(Pc(1,j-1)-Pc(1,j))))...
               -Po(1,j)*(Cpoo(1,j)+aalpha(1,j)*Cpog(1,j))+Qsco(SSS)+aalpha(1,j)*Qsco(SSS)*((Landag(1,j)/Landao(1,j))+Rs(1,j));
               end
           

               for j=2:gte-1
               F(gr+1,j)=-(aalpha(gr+1,j)*(TrgMinus(gr+1,j)*(Pc(gr,j)-Pc(gr+1,j))...
               +TtegPlus(gr+1,j)*(Pc(gr+1,j+1)-Pc(gr+1,j))...
               +TtegMinus(gr+1,j)*(Pc(gr+1,j-1)-Pc(gr+1,j))))...
               -Po(gr+1,j)*(Cpoo(gr+1,j)+aalpha(gr+1,j)*Cpog(gr+1,j));
               end
               
           for i=2:gr
           F(i,1)=-(aalpha(i,1)*(TrgPlus(i,1)*(Pc(i+1,1)-Pc(i,1))+...
               TrgMinus(i,1)*(Pc(i-1,1)-Pc(i,1))+TtegPlus(i,1)*(Pc(i,2)-Pc(i,1))...
               +TtegMinus(i,1)*(Pc(i,gte)-Pc(i,1))))...
               -Po(i,1)*(Cpoo(i,1)+aalpha(i,1)*Cpog(i,1));
           end
           
           
           for i=2:gr
               F(i,gte)=-(aalpha(i,gte)*(TrgPlus(i,gte)*(Pc(i+1,gte)-Pc(i,gte))+...
               TrgMinus(i,gte)*(Pc(i-1,gte)-Pc(i,gte))+TtegPlus(i,gte)*(Pc(i,1)-Pc(i,gte))...
               +TtegMinus(i,gte)*(Pc(i,gte-1)-Pc(i,gte))))...
               -Po(i,gte)*(Cpoo(i,gte)+aalpha(i,gte)*Cpog(i,gte));
           end
  
           F(gr+1,1)=-(aalpha(gr+1,1)*(TrgMinus(gr+1,1)*(Pc(gr,1)-Pc(gr+1,1))...
               +TtegPlus(gr+1,1)*(Pc(gr+1,2)-Pc(gr+1,1))...
               +TtegMinus(gr+1,1)*(Pc(gr+1,gte)-Pc(gr+1,1))))...
               -Po(gr+1,1)*(Cpoo(gr+1,1)+aalpha(gr+1,1)*Cpog(gr+1,1));

           F(gr+1,gte)=-(aalpha(gr+1,gte)*(TrgMinus(gr+1,gte)*(Pc(gr,gte)-Pc(gr+1,gte))...
               +TtegPlus(gr+1,gte)*(Pc(gr+1,1)-Pc(gr+1,gte))...
               +TtegMinus(gr+1,gte)*(Pc(gr+1,gte-1)-Pc(gr+1,gte))))...
               -Po(gr+1,gte)*(Cpoo(gr+1,gte)+aalpha(gr+1,gte)*Cpog(gr+1,gte));
              
           F(1,1)=-aalpha(1,1)*((TrgPlus(1,1)*(Pc(2,1)-Pc(1,1))+...
               TtegPlus(1,1)*(Pc(1,2)-Pc(1,1))...
               +TtegMinus(1,1)*(Pc(1,gte)-Pc(1,1))))...
               -Po(1,1)*(Cpoo(1,1)+aalpha(1,1)*Cpog(1,1))+Qsco(SSS)+aalpha(1,1)*Qsco(SSS)*((Landag(1,1)/Landao(1,1))+Rs(1,1));

 
           F(1,gte)=-(aalpha(1,gte)*(TrgPlus(1,gte)*(Pc(2,gte)-Pc(1,gte))+...
               TtegPlus(1,gte)*(Pc(1,1)-Pc(1,gte))...
               +TtegMinus(1,gte)*(Pc(1,gte-1)-Pc(1,gte))))...
               -Po(1,gte)*(Cpoo(1,gte)+aalpha(1,gte)*Cpog(1,gte))+Qsco(SSS)+aalpha(1,gte)*Qsco(SSS)*((Landag(1,gte)/Landao(1,gte))+Rs(1,gte));   
           
           CC=zeros((gr+1)*gte);% Matrix Zarayeb
           for i=1:(gr+1)*gte
               CC(i,i)=B(i);
               if  i<=((gr+1)*gte-1)
                   CC(i,i+1)=C(i);
                   CC(i+1,i)=A(i+1);
               end
               if i>=1 && i<=((gr+1)*gte-(gr+1))
                   CC(i,i+(gr+1))=D(i);
                   CC(i+(gr+1),i)=E(i+(gr+1));
               end
               if i>((gr+1)*gte-(gr+1)) 
                   CC(i,i-((gr+1)*gte-(gr+1)))=D(i);
                   CC(i-((gr+1)*gte-(gr+1)),i)=E(i-((gr+1)*gte-(gr+1)));
               end
           end
           
          
      
           for i=1:(gr+1)*gte
               BB(i)=F(i);
           end

           
           BBB=BB';
           PM=CC\BBB;
           

           
           for i=1:gte
               Po(:,i)=PM(1+(i-1)*(gr+1):i*(gr+1));
           end
           

           
%--------------------------gas Saturation-------------------------------%
%=========================================================================%       
 
           for i=2:gr
               for j=2:gte-1
           Sg(i,j)=Sgg(i,j)+(1/Csgo(i,j))*(TroPlus(i,j)*(Po(i+1,j)-Po(i,j))+TroMinus(i,j)*(Po(i-1,j)-Po(i,j))+TteoPlus(i,j)*(Po(i,j+1)-Po(i,j))+TteoMinus(i,j)*(Po(i,j-1)-Po(i,j))-Cpoo(i,j)*(Po(i,j)-Poo(i,j)));
               end
           end
           
           for j=2:gte-1
           Sg(1,j)=Sgg(1,j)+(1/Csgo(1,j))*(TroPlus(1,j)*(Po(2,j)-Po(1,j))+TteoPlus(1,j)*(Po(1,j+1)-Po(1,j))+TteoMinus(1,j)*(Po(1,j-1)-Po(1,j))-Qsco(SSS)-Cpoo(1,j)*(Po(1,j)-Poo(1,j)));
           end
           
          
           for j=2:gte-1
           Sg(gr+1,j)=Sgg(gr+1,j)+(1/Csgo(gr+1,j))*(TroMinus(gr+1,j)*(Po(gr,j)-Po(gr+1,j))+TteoPlus(gr+1,j)*(Po(gr+1,j+1)-Po(gr+1,j))+TteoMinus(gr+1,j)*(Po(gr+1,j-1)-Po(gr+1,j))-Cpoo(gr+1,j)*(Po(gr+1,j)-Poo(gr+1,j)));
           end
          
           for i=2:gr 
           Sg(i,1)=Sgg(i,1)+(1/Csgo(i,1))*(TroPlus(i,1)*(Po(i+1,1)-Po(i,1))+TroMinus(i,1)*(Po(i-1,1)-Po(i,1))+TteoPlus(i,1)*(Po(i,2)-Po(i,1))+TteoMinus(i,1)*(Po(i,gte)-Po(i,1))-Cpoo(i,1)*(Po(i,1)-Poo(i,1)));
           end
           
           for i=2:gr
           Sg(i,gte)=Sgg(i,gte)+(1/Csgo(i,gte))*(TroPlus(i,gte)*(Po(i+1,gte)-Po(i,gte))+TroMinus(i,gte)*(Po(i-1,gte)-Po(i,gte))+TteoPlus(i,gte)*(Po(i,1)-Po(i,gte))+TteoMinus(i,gte)*(Po(i,gte-1)-Po(i,gte))-Cpoo(i,gte)*(Po(i,gte)-Poo(i,gte)));
           end
           
           Sg(gr+1,1)=Sgg(gr+1,1)+(1/Csgo(gr+1,1))*(TroMinus(gr+1,1)*(Po(gr,1)-Po(gr+1,1))+TteoPlus(gr+1,1)*(Po(gr+1,2)-Po(gr+1,1))+TteoMinus(gr+1,1)*(Po(gr+1,gte)-Po(gr+1,1))-Cpoo(gr+1,1)*(Po(gr+1,1)-Poo(gr+1,1)));
           
           Sg(gr+1,gte)=Sgg(gr+1,gte)+(1/Csgo(gr+1,gte))*(TroMinus(gr+1,gte)*(Po(gr,gte)-Po(gr+1,gte))+TteoPlus(gr+1,gte)*(Po(gr+1,1)-Po(gr+1,gte))+TteoMinus(gr+1,gte)*(Po(gr+1,gte-1)-Po(gr+1,gte))-Cpoo(gr+1,gte)*(Po(gr+1,gte)-Poo(gr+1,gte)));
           
           Sg(1,1)=Sgg(1,1)+(1/Csgo(1,1))*(TroPlus(1,1)*(Po(2,1)-Po(1,1))+TteoPlus(1,1)*(Po(1,2)-Po(1,1))+TteoMinus(1,1)*(Po(1,gte)-Po(1,1))-Qsco(SSS)-Cpoo(1,1)*(Po(1,1)-Poo(1,1)));
           
           Sg(1,gte)=Sgg(1,gte)+(1/Csgo(1,gte))*(TroPlus(1,gte)*(Po(2,gte)-Po(1,gte))+TteoPlus(1,gte)*(Po(1,1)-Po(1,gte))+TteoMinus(1,gte)*(Po(1,gte-1)-Po(1,gte))-Qsco(SSS)-Cpoo(1,gte)*(Po(1,gte)-Poo(1,gte)));



           So=1-Sg-Swi;
           

for xx=1:gr+1
    for yy=1:gte
        if Sg(xx,yy)>=(1-Sor-Swi)
            Sg(xx,yy)=1-Sor-Swi;
        end
    end
end

SSg(iii,SSS)=Sg(1,gte);

for jj=1:gte
PwfPro(iii,SSS)=Po(2,jj)-Qsco(SSS)*(Muo(2,jj)*Bo(2,jj)*log(Re/Rw))/(7.08*K*Kro(2,jj)*deltaz);
end


end
SSg(2:NT(SSS)+1,SSS)=SSg(1:NT(SSS),SSS);
PwfPro(2:NT(SSS)+1,SSS)=PwfPro(1:NT(SSS),SSS);
PwfPro(1,1)=Poi;



for i=2:SSS
SSg(1,SSS)=SSg(NT(SSS-1)+1,SSS-1);   
PwfPro(1,SSS)=PwfPro(NT(SSS-1)+1,SSS-1);
end


end
%--------------------------------Pwf Plot---------------------------------%
%=========================================================================%
SSg
PwfPro
for i=2:length(NT)+1
    SumNT(i)=sum(NT(1:i-1));
end
SumNT(1)=0;

for j=1:SSS
    subplot(2,1,1)
    plot(SumNT(j):SumNT(j+1),PwfPro(1:NT(j)+1,j),'r')
    hold on
    xlabel('Time (hours)')
    ylabel('Pwf-Production Well')
    grid minor
    grid on
    subplot(2,1,2)
    plot(linspace(SumNT(j),SumNT(j+1),10000),Qsco(j),'*b')
    hold on
    xlabel('Time (hours)')
    ylabel('Oil Rate')
    grid minor
    grid on
    
end
figure;
for j=1:SSS

    plot(SumNT(j):SumNT(j+1),SSg(1:NT(j)+1,j),'r')
    hold on
    xlabel('Time (hours)')
    ylabel('Gas Saturation(%)')
    grid minor
    grid on  
end
%-------------------------------ERROR-------------------------------------%
%=========================================================================%
IOIP=Phi*deltaz*(Re^2)*(1-Swi)/(5.615*(0.9759+0.00012*((Gamag*((Poi/18.2+1.4)*10^(0.0125*API-0.0009*(TR-460)))^1.2048)*(Gamag/Gamao)^0.5+1.25*(TR-460))^1.2));
for i=1:gr+1
 NN(i)=Phi.*Vb(i).*So(i)./(5.615*Bo(i));
end
Nr=sum(NN);
AAAAA=IOIP-Nr;
disp(['(IOIP - Residual Oil) = ',num2str(AAAAA)])
for i=1:TTD
DEBI(i)=Qsco(i)*NT(i)/24;
end
BBBB=sum(DEBI);
disp(['Sum of Oil Rates = ',num2str(BBBB)])
ER=abs(AAAAA-BBBB)/AAAAA;
disp(['Error = ',num2str(ER)])
