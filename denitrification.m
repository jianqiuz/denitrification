function dcdt=denitrification(t,y)
global par

k_NO3=0.013;%mol/L% reaction rate constant of nitrate reduction ##d[NO3]/dt=-2*r_NO3
k_NO2=0.00088;%mol/L% reaction rate constant of nitrite reduction##d[NO2]/dt=2*r_NO3-4*r_NO2
k_NO=8.5*10^(-9);%mol/L% reaction rate constant of NO reduction##d[NO]/dt=4*r_NO2-8*r_NO
k_N2O=5*10^(-6);%mol/L% reaction rate constant of N2O reduction##d[N2O]/dt=4*r_NO-4*r_N2O, d[N2]/dt=4*r_N2O
%%Michaelis-Menten constants from Gu's 
%%%All microbially mediated transformations of N species are aqueous
%%%reactions, denoted by C_NO3,C_NO2,C_NO,C_N2O,C_N2.
k_OC=1.42*10^(-3); %mol/L, from (Li,2000) half-saturation constant, used for Monod kinetics of microbial biomass growing
kI_NO=0.00002; %mol/L, inhibition coefficient by NO, Andrew-Haldane model applied ( growth kinetics of inhibitory substrate)
k_O2=2.52*10^(-5);
k_I_NO_O2=0.0174*10^(-6);

Vm_NO3=0; %h-1 maximum enzyme synthesis rate
Vm_NO2=0; %h-1 mimum enzyme synthesis rate
Vm_NO=0; %h-1 maximum enzyme synthesis rate
Vm_N2O=0; % smaller enzyme synthesis rate for NOS
 
Ke_NO3=10^(-11); % mol/L  enzyme saturation coefficient
Ke_NO2=5*10^(-5);
Ke_NO=5.4*10^(-8);
Ke_N2O=5*10^(-7);%%enzyme saturation coefficient
 
KIO_NO3=2.5*10^(-5);
KIO_NO2=2.2*10^(-5);
KIO_NO=4*10^(-4);
KIO_N2O=1*10^(-7);



Yc=0.503; %biomass yeild coefficient( growing upon carbon) ~~~carbon utilization rate, from Li,1992


%%Henry's law constant
kla=19.3; %h-1 mass transfer coefficient
%H_CO2=1.2;
H_N2O=1.64;%%henry's law constant gas/liquid dementionless (0.025 M/atm)
H_NO=21.5;%%henry's law constant gas/liquid dementionless (0.0019 M/atm)
H_N2=64.9;%%henry's law constant gas/liquid dementionless (0.00063 M/atm)
H_O2=33;

V_gas=0.035; %L
V_liquid=0.005; %L
Vr=V_liquid/V_gas;

%initial net reaction rates for each state varible
r_NO3=0; %mol/L/h
r_NO2=0;
r_NO=0;
r_N2O=0;
r_N2=0;
r_OC=0;
r_O2=0;
r_B=0;
r_E_nar=0;
r_E_nir=0;
r_E_nor=0;
r_E_nos=0;


C_NO3=y(1);%%%aqueous phase concentration in mol/L 
C_NO2=y(2);
C_NO=y(3);
C_N2O=y(4);
C_N2=y(5);
C_OC=y(6);
C_NOg=y(7);
C_N2Og=y(8);
C_N2g=y(9);
C_O2=y(10);
C_O2g=y(11);
C_B=y(12);
E_nar=y(13);
E_nir=y(14);
E_nor=y(15);
E_nos=y(16);


%%reactions%%
%%(max growth rate on the substrates)*(reaction dynamic)*(Biomass)*(Enzymes)


g_O2=0.2*C_B*(C_O2/(k_O2*(1+C_NO/k_I_NO_O2)+C_O2))*(C_OC/(k_OC+C_OC));
%g_O2=0;
r_O2=r_O2-g_O2;

%g_NO3=0;
g_NO3=0.648*C_B*0.2*(C_NO3/(k_NO3+C_NO3))*E_nar*(C_OC/(k_OC+C_OC));%%growing reaction rate*(0.000004/(C_O2+0.000004))*(C_O2/(C_O2+0.000004))*(C_O2/(C_O2+0.000004))
r_NO3=r_NO3-2*g_NO3;
r_NO2=r_NO2+2*g_NO3;


g_NO2=0.648*C_B*0.2*(C_NO2/(C_NO2+k_NO2))*E_nir*(C_OC/(k_OC+C_OC));
r_NO2=r_NO2-4*g_NO2;
r_NO=r_NO+4*g_NO2;


g_NO=0.3265*C_B*0.2*(C_NO/(C_NO*(1+C_NO/kI_NO)+k_NO))^2*E_nor*(C_OC/(k_OC+C_OC)); %(C_NO/(k_NO+C_NO))^2   *(C_OC/(k_OC+C_OC))
r_NO=r_NO-8*g_NO;
r_N2O=r_N2O+4*g_NO;

g_N2O=0.3247*C_B*0.2*(C_N2O/(k_N2O+C_N2O))*E_nos*(C_OC/(k_OC+C_OC));
%g_N2O=0;
r_N2O=r_N2O-4*g_N2O;
r_N2=r_N2+4*g_N2O;


%r_B=0;
r_B=(r_B+g_O2+g_NO3+g_NO2+2*g_NO+2*g_N2O)*0.3*5-0.001*C_B;%%%problem CO2 production or biomass incoporation, use Yc to illustr
%%m%%aximum growth rate on different substrate:0.063/0.063/0.02988/0.030132

r_OC=r_OC-g_O2-g_NO2-2*g_NO-2*g_N2O;
%r_OC=r_OC-g_O2-3*r_B;
    


r_E_nar=Vm_NO3*(C_NO3/(Ke_NO3+C_NO3))*(1-E_nar)*(KIO_NO3/(C_O2+KIO_NO3));
r_E_nir=Vm_NO2*(C_NO2/(Ke_NO2+C_NO2))*(1-E_nir)*(KIO_NO2/(C_O2+KIO_NO2));
r_E_nor=Vm_NO*(C_NO/(Ke_NO+C_NO))*(1-E_nor)*(KIO_NO/(C_O2+KIO_NO));
r_E_nos=Vm_N2O*(C_N2O/(Ke_N2O+C_N2O))*(1-E_nos)*(KIO_N2O/(C_O2+KIO_N2O));




%%gas/liquid transfer rates
rtr_NO=kla*(C_NOg/H_NO-C_NO);
rtr_N2O=kla*(C_N2Og/H_N2O-C_N2O);
rtr_N2=kla*(C_N2g/H_N2-C_N2);
rtr_O2=kla*(C_O2g/H_O2-C_O2);


%Initialize vector of detrivatives
dcdt=zeros(16,1);
dcdt(1)=r_NO3;
dcdt(2)=r_NO2;
dcdt(3)=r_NO+rtr_NO;
dcdt(4)=r_N2O+rtr_N2O;
dcdt(5)=r_N2+rtr_N2;
dcdt(6)=r_OC;%%organic carbon
dcdt(7)=-Vr*rtr_NO; %NO gas
dcdt(8)=-Vr*rtr_N2O; %N2O gas
dcdt(9)=-Vr*rtr_N2;
dcdt(10)=r_O2+rtr_O2;
dcdt(11)=-Vr*rtr_O2;
dcdt(12)=r_B;
dcdt(13)=r_E_nar;
dcdt(14)=r_E_nir;
dcdt(15)=r_E_nor;
dcdt(16)=r_E_nos;


end


