function [DRAG] = Drag(WING,V_TAIL,FUSELAGE,PARAM)

M = PARAM.M; % mach number
S = WING.Sw; % Surface wing ;
Snet = WING.Sw - FUSELAGE.a_el*2 * WING.cw_root; % Surface net (wing-fusellage) ??????????
CL0 = 0.13;% CL at 0 degree angle of attack 
CLa = WING.CLw_alpha; % CL alpha 
A = WING.A; %AR of the wing
et = abs(WING.theta_tip *180/pi); %Twist angle [deg](doit ??tre positif)
tc = 0.12; %Thickness to chord ratio (depend du profil?? d'airfoil choisi)
Sh = V_TAIL.Sh_tail; %Horizontal surface of the tail 
Ah = V_TAIL.AR_T; %Aspect ratio of the horizontal tail
CLh = V_TAIL.CL_tail; % Tail plane lift coeff
tch = 0.1; % Htail thickness ratio
sweep12h = V_TAIL.Lambda_T; %Htail mid sweep angle [deg]
sweeph = V_TAIL.Lambda_T; %Htail leading edge sweep angle [deg}
tcv = 0.1; %Vtail thick ratio
sweep12v = V_TAIL.Lambda_T; % Vtail mid sweep angle
Sv = V_TAIL.Sv_tail; % Vertical tail surface
cr = WING.cw_root; % root chord wing 
Vf = 8.5; % Volume of the fuselage ->premi??re approx voir CAD, prev = 6
Ac = FUSELAGE.a_el * FUSELAGE.b_el * pi; %Cross section area 
lf = FUSELAGE.l_f; %fuselage length 
Sfwet = 30; %fuselage wetted area -> CAD, prev = 24
Snwet = 4.6/2; %nacelle wetted area -> CAD, prev = 5.5/2
Spwet = 0.8*pi*1.57; %engine wetted area (disq of engine) ->>>>>>>>>>> ENGINE CARACT
Mg = 0.7; %fully extended flow mach number p507
Minf = 0.7; %Flying mach number
nf = 2*FUSELAGE.a_el/WING.bw; %fuselage diameter/wingspan
taper = WING.tap; %wing taper ratio 
sweep12 = WING.sweep * 180/pi; %wing mid sweep angle
tr = 0.12 * WING.cw_root; % root thickness p509
sweep = WING.sweep; %leading edge sweep angle [deg]
deda = PARAM.de_dAOA; %Downash slope --<><<<<<<<<<<<<<<<<<<<<<>>>>> ligne 287 AIRCRAFT
l = WING.cw_MAC; %length p500 (MAC je suppose) 
SMC = WING.Sw /WING.bw; %Standard mean chord
airEntryLength = 1.6;  %Longeur entre d'air -> CAD, prev = 3.98
ch = V_TAIL.hMAC; %Htail mean chord ->rajouter depuis aircraft ligne 119
cv = V_TAIL.vMAC; %Vtail mean chord ->rajouter depuis aircraft ligne 128
e = 0.8; % Coef d'oswald ?????? 
CLi = 0.409; % J'ai pas bien compris (mis ?? 0.36 dnas le code de l??o)
clmax = 1.4; %PARAM.CL; % DEMANDER LEO 
CLmax = 1.65; %PARAM.CL; % DEMANDER LEO 
CL = PARAM.CL; %voir si c'est l? le probl?me
beta = sqrt(1-Minf^2);
CD0 = WING.CD0;

%From graph 
delta = 0.0014;% Voir figure 1 page 514 du pdf
Col = ((1.4*10^(-5) + 1.15*10^(-5))/(2*beta));%Voir figure 2 page 515 du pdf
Cll = -0.25*10^(-4)/beta; % Voir figure 3 page 515 du pdf
CDB = 0.02;%Paremetre voir graphique au dessus graphique et FR se calcule avec les dimensions du moteur
B = 0; %Figure 13 page 524 du pdf ????????????????????????
A1 = 4.8; %Figure 13 page 524 du pdf ??????????????????????, prev = 4
A2 = 5.26; %Figure 13 page 524 du pdf ???? prev = 4.55
Cci = 4 *cr; %Page 510 du livre ???????????????????????
lN = 3; %Figure 10 page 524 du PDF + CAD, prev = 2.36
lA = 4.191; %Figure 10 page 524 du PDF + CAD 

%% Vortex
CDvF9 = (1 + delta)*(CL^2)/pi/A;
deCdvF21 = et^2 * Col + et * CL * Cll;
af = (CL - CL0)/CLa; %in rad
CDSF24 = 0.15 * af^2 * Vf^(2/3);
CDSF26 = 1.02 * CLh^2 * Sh / (pi*Ah);
CDVortexWing = CDvF9 + deCdvF21;
CDVortexFus = CDSF24/S;
CDVortexTail = CDSF26/S;
CDVortex = CDVortexWing + CDVortexFus +CDVortexTail;

%% Profile drag
phiw = 2.7 * tc + 100 * tc^4;
cdpminF35 = 2 * getCF(real(SMC),M) * (1 + phiw *cosd(sweep12)^2);
dlcdpref = 0.1 * clmax - 0.0046 * (l + 2.75 *tc);
CDpF36 = cdpminF35 * Snet/S + 0.75 * dlcdpref * ((CL-CLi)/(CLmax-CLi))^2;
Dfeff = sqrt(4/pi*Ac);
sigmaeff = min(lf/Dfeff,(lN+lA)/Dfeff + 2);
phif = 2.2/(sigmaeff^1.5) + 3.8/(sigmaeff^3);% bcse lA/Dfeff > 2 (see p.523 pdf)
CDSF = getCF(real(lf),M) * Sfwet * (1 + phif);
CDSbasicF46 = CDSF + CDB * pi/4 * Dfeff^2;
dabCDSF51 = A1 * abs(sin(af)^3) + A2 * abs(sin(af-B)^3)/cos(B);
CDSF58 = 1.25 * getCF(real(airEntryLength),M) * Snwet;
CDSF60 = getCF(real(airEntryLength),M)*Spwet*(Mg/Minf)^(11/6) * ((1+0.116*Minf^2)/(1+0.116*Mg^2))^(2/3);
CDShbasicF63 = 2*getCF(real(ch),M)*(1+2.75*tch*cosd(sweep12h)^2)*Sh;
dlCDShF64 = 0.33*CLh^2 * Sh / (cosd(sweeph)^2*pi*Ah);
CDSvF56 = 2*getCF(real(cv),M)*(1+2.75*tcv*cosd(sweep12v)^2)*Sv;

CDProfileWing = CDpF36;
CDProfileFuselage = (CDSbasicF46+dabCDSF51)/S;
CDProfileNacelles = CDSF58/S + CDSF60/S;
CDProfileNacelles = CDProfileNacelles*2;
CDProfileHTail = (CDShbasicF63 + dlCDShF64)/S;
CDProfileVTail = CDSvF56/S;
CDProfileTail = CDProfileHTail+ CDProfileVTail;
CDProfileTotal = CDProfileWing + CDProfileFuselage + CDProfileNacelles + CDProfileHTail + CDProfileVTail;

%% Interference corrections
diCDvF65 = 0.55*nf/(1+taper)*(2-pi*nf)*CL0^2/(pi*A);
diCDSpF66 = 1.5*getCF(SMC,M)*tr*Cci*cosd(sweep12)^2;
%afp = B*(sqrt(A1/A2)-1)/(A1/A2-1);
afp = 0;
D1 = A1 * abs(sin(afp)^3) + A2 * abs(sin(afp-B)^3)/cos(B);
diCDSF67 = -0.81*getCF(SMC,M)*CL*cr*D1;
dedCL = deda/CLa; %je crois
diCDSvF71 = CLh*CL*(dedCL - 2/(pi*A))*Sh;

CDInterference = (diCDvF65 + diCDSpF66 + diCDSF67 + diCDSvF71)/S;

%% Protuberance, surface impefection
CdImpWing = 0.06*CDProfileWing;
CdImpFusEmp = 0.07*(CDProfileFuselage + CDProfileHTail + CDProfileVTail);
CdImpSyst = 0.03*CD0;

CDImperfection = CdImpWing + CdImpFusEmp + CdImpSyst;
CDOther = CDInterference + CDImperfection;

CDTotal = CDVortex + CDProfileTotal + CDInterference + CDImperfection;

DRAG = table(CDTotal,CDProfileTotal,CDVortex,CDOther,CDProfileWing,CDProfileFuselage,CDProfileNacelles,CDProfileTail);

end 