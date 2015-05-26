function GillespieToyNetwork

close all;

uncatalysed = 0.1;
catalysed = 1.0;

k1 = uncatalysed;
k2 = catalysed;
k3 = uncatalysed;
k4 = catalysed;
k5 = uncatalysed;
k6 = uncatalysed;
k7 = uncatalysed;
k8 = uncatalysed;
k9 = uncatalysed;
k10 = uncatalysed;
k11 = uncatalysed;
k12 = uncatalysed;
k13 = uncatalysed;
k14 = uncatalysed;

k_values = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14];

O_initial = 10;
CN_initial = 10;
C3N_initial = 10;
CO_initial = 10;
N_initial = 10;
NC2N_initial = 0;
NC4N_initial = 0;
OCCO_initial = 0;
NCCO2_initial = 0;

%stoichiometric_matrix = NEED TO ADD STOICHIOMETRIC VALUES... BE CAREFUL
%WHEN TYPING! MOLECULES SHOULD BE ROWS, REACTION SHOULD BE COLUMNS.

realisations = 10;

for i = 1:realisations

    O = O_initial;
    N = N_initial;
    CO = CO_initial;
    CN = CN_initial;
    C3N = C3N_initial;
    NC2N = NC2N_initial;
    NC4N = NC4N_initial;
    OCCO = OCCO_initial;
    NCCO2 = NCCO2_initial;
    
    molecule_numbers = [O N CO CN C3N NC2N NC4N OCCO NCCO2];
    
    time=0;
    k=1;
    
    Oplot(1,i)=O;
    Nplot(1,i)=N;
    COplot(1,i) = CO;
    CNplot(1,i)=CN;
    C3Nplot(1,i) = C3N;
    NC2Nplot(1,i)=NC2N;
    NC4Nplot(1,i) = NC4N;
    OCCOplot(1,i) = OCCO;
    NCCO2plot(1,i) = NCCO2;
    
    timeplot(1,i)=0;
    
    while k < 10000
        
        rr = rand(2,1);
        
        prop(1) = CN * (CN - 1) * k_values(1);
        prop(2) = CN * C3N * (k_values(2) * C3N * NC2N); %This reaction is catalysed by C3N and NC2N. kinetic constant is multiplied by catalyst numbers, as in Hordijk et al. (2014).
        prop(3) = CN * O * CO * N * k_values(3);
        prop(4) = CO * (CO - 1) * (k_values(4) * NC4N); % Catalysed by NC4N.
        prop(5) = k_values(5);
        prop(6) = k_values(6);
        prop(7) = k_values(7);
        prop(8) = k_values(8);
        prop(9) = k_values(9);
        prop(10) = NC4N * k_values(10);
        prop(11) = NC2N * k_values(11);
        prop(12) = N * k_values(12);
        prop(13) = NCCO2 * k_values(13);
        prop(14) = OCCO * k_values(14);
        
        alpha0 = sum(prop);
        tau = (1/alpha0)*log(1/rr(1));
        time = time + tau;
        
        CumulativePropSum = cumsum(prop);
        
        for i = 1:length(prop)
            if rr(2) * alpha0 >= CumulativePropSum(i) && i < length(prop)
                if rr(2) * alpha0 < CumulativePropSum(i+1)
                    %CARRY OUT REACTION i+1...
                    
                
        
        
        
end
