%%%denitrification modeling%%%
%%NO3/NO2/NO/N2O/N2
global k_NO3 k_NO2 k_NO k_N2O k_OC kI_NO k_O2 k_I_NO_O2

k_NO3=0.013;%mol/L% reaction rate constant of nitrate reduction ##d[NO3]/dt=-2*r_NO3
k_NO2=0.00088;%mol/L% reaction rate constant of nitrite reduction##d[NO2]/dt=2*r_NO3-4*r_NO2
k_NO=8.5*10^(-9);%mol/L% revised according to conrad
k_N2O=5*10^(-6);%mol/L% from conrad
%%Michaelis-Menten constants from Gu's 
%%%All microbially mediated transformations of N species are aqueous
%%%reactions, denoted by C_NO3,C_NO2,C_NO,C_N2O,C_N2.
k_OC=1.42*10^(-3); %mol/L, from Maggi, gu (2009)
kI_NO=0.00002; %mol/L, inhibition coefficient by NO, Andrew-Haldane model applied growth kinetics of inhibitory substrate)
k_O2=2.52*10^(-5);%check
k_I_NO_O2=0.0174*10^(-6);

%%initial condition (concentration at time 0)
C0_NO3=0.00005; %%mol/L
C0_NO2=0;
C0_NO=0;
C0_N2O=0;
C0_N2=0;
C0_OC=0.0185;%%mol/L
C0_NOg=0;
C0_N2Og=0;
C0_N2g=0;
C0_O2=0.0000124;%mol/L aqueous phase
C0_O2g=0.000040;%mol/L gas phase
C0_B=0.00314; %initial active biomass 291 umol PLFAs /L, from PLFA measurements
E0_nar=0.42;
E0_nir=0.42;
E0_nor=0.42;
E0_nos=0;



C0=[C0_NO3 C0_NO2  C0_NO  C0_N2O C0_N2 C0_OC C0_NOg C0_N2Og C0_N2g C0_O2 C0_O2g C0_B E0_nar E0_nir E0_nor E0_nos];% p_N2O r_N2O];%%% E_nar  E_nir E_nor E_nos C_B_NO3 C_B_NO2 C_B_NO C_B_N2O];


options=odeset('RelTol',[1e-9],'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]);
[T,Y]=ode23(@denitrification,[0 48],C0,options);

%plot(T,Y(:,7),'s',T,Y(:,8),'d',T,Y(:,9),'x');
%plot(T,Y(:,13),'s',T,Y(:,14),'d',T,Y(:,15),'x',T,Y(:,16),'o');
%plot(T,Y(:,1),'-',T,Y(:,2),'o',T,Y(:,3),'s',T,Y(:,4),'d',T,Y(:,5),'x');
%plot(T,Y(:,2),'o',T,Y(:,3),'s',T,Y(:,4),'d',T,Y(:,5),'x');

%plot(T,Y(:,13),'s',T,Y(:,14),'d',T,Y(:,15),'x',T,Y(:,16),'o');
%plot(T,Y(:,10),'x',T,Y(:,11),'o');
%plot(T,Y(:,4),'x',T,Y(:,8),'o');
%plot(T,Y(:,12));


H=1000000000*Y; %%Y(mol/L) to nmol/L


%soil incubation 0-5 cm

%plot(T,H(:,8),0,1,'s',3,150,'s',6,219,'s',12,629,'s',24,1434,'s',36,1952,'s',48,2128,'s');%CHL+C2H2
%plot(T,H(:,8),0,1,'s',3,30,'s',6,140,'s',12,567,'s',24,1275,'s',36,1862,'s',48,2022,'s');%CHL

%plot(T,H(:,8),0,1.6,'s',3,298,'s',6,696,'s',12,1497,'s',24,2661,'s',36,2837,'s',48,2758,'s');%C2H2

%plot(T,H(:,8),0,1.3,'s',3,79,'s',6,526,'s',12,1046,'s',24,273,'s',36,19,'s',48,9,'s');%H2O'

%plot(T,H(:,9),0,0.4,'s',3,215,'s',6,170,'s',12,451,'s',24,2334,'s',36,2816,'s',48,2747,'s');%N2
%production 


%%%soil incubation 5-10 cm

%plot(T,H(:,8),0,1,'s',3,55,'s',6,132,'s',12,210,'s',24,397,'s',36,474,'s',48,637,'s');%CHL+C2H2

%plot(T,H(:,8),0,2.24,'s',3,100,'s',6,223,'s',12,439,'s',24,753,'s',36,1049,'s',48,1366,'s');%C2H2

%plot(T,H(:,8),0,1.6,'s',3,58,'s',6,121,'s',12,267,'s',24,84,'s',36,15,'s',48,9,'s');%H2O

%plot(T,H(:,9),0,0.7,'s',3,43,'s',6,102,'s',12,172,'s',24,669,'s',36,1034,'s',48,1357,'s');%N2



%%%soil incubation 10-15 cm

%plot(T,H(:,8),0,2,'s',3,5,'s',6,16,'s',12,87,'s',24,200,'s',36,318,'s',48,399,'s');%CHL+C2H2

%plot(T,H(:,8),0,3,'s',3,5,'s',6,82,'s',12,412,'s',24,686,'s',36,893,'s',48,1085,'s');%C2H2

%plot(T,H(:,8),0,1.2,'s',3,4,'s',6,39,'s',12,275,'s',24,217,'s',36,120,'s',48,71,'s');%H2O

%plot(T,H(:,9),0,1,'s',3,1.6,'s',6,43,'s',12,158,'s',24,469,'s',36,773,'s',48,1014,'s');%N2

%plot(T,H(:,8),0,1,'s',3,62,'s',6,105,'s',12,157,'s',24,136,'s',36,126,'s',48,66,'s');%CHL

%%%soil incubation 15-25 cm

%plot(T,H(:,8),0,2.7,'s',3,5.4,'s',6,45,'s',12,154,'s',24,312,'s',36,445,'s',48,579,'s');%CHL+C2H2

%plot(T,H(:,8),0,2,'s',3,7,'s',6,45,'s',12,252,'s',24,650,'s',36,858,'s',48,952,'s');%C2H2

%plot(T,H(:,8),0,1,'s',3,5,'s',6,42,'s',12,244,'s',24,501,'s',36,375,'s',48,252,'s');%H2O

%plot(T,H(:,9),0,0.5,'s',3,2,'s',6,5,'s',12,9,'s',24,148,'s',36,483,'s',48,701,'s');%N2

%plot(T,H(:,8),0,2,'s',3,3.5,'s',6,18,'s',12,53,'s',24,120,'s',36,183,'s',48,198,'s');%CHL


