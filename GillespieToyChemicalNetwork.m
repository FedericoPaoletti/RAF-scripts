function GillespieToyNetwork 

close all;
%rand('state',8);

uncatalysed = 0.1; %
catalysed = 1.0;
in = 0.3;
out = 0.01;

k1 = uncatalysed;
k2 = uncatalysed;
k3 = uncatalysed;
k4 = uncatalysed;
k5 = in;
k6 = in;
k7 = in;
k8 = in;
k9 = in;
k10 = out;
k11 = out;
k12 = out;
k13 = out;
k14 = out;
%k15 = out;

k_values = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14];

O_initial = 10;
N_initial = 10;
CO_initial = 10;
CN_initial = 10;
C3N_initial = 10;
NC2N_initial = 0;
NC4N_initial = 0;
OCCO_initial = 0;
NCCO2_initial = 0;

stoichiometric_matrix = [0 0 -1 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 -1 0 0; 0 0 -1 -2 0 0 1 0 0 0 0 0 0 0; -2 -1 -1 0 1 0 0 0 0 0 0 0 0 0; 0 -1 0 0 0 0 0 0 1 0 0 0 0 0; 1 0 0 0 0 0 0 0 0 0 -1 0 0 0; 0 1 0 0 0 0 0 0 0 -1 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 -1; 0 0 1 0 0 0 0 0 0 0 0 0 -1 0];
stoichiometric_matrix;                            

realisations = 1;

for i = 1:realisations

    molecules = [O_initial N_initial CO_initial CN_initial C3N_initial NC2N_initial NC4N_initial OCCO_initial NCCO2_initial];
    
    time = 0;
    k = 1;
    
    Oplot(1,i) = O_initial;
    Nplot(1,i) = N_initial;
    COplot(1,i) = CO_initial;
    CNplot(1,i)= CN_initial;
    C3Nplot(1,i) = C3N_initial;
    NC2Nplot(1,i)= NC2N_initial;
    NC4Nplot(1,i) = NC4N_initial;
    OCCOplot(1,i) = OCCO_initial;
    NCCO2plot(1,i) = NCCO2_initial;
    
    timeplot(1,i) = 0;
    
    while k < 10000
        
        rr = rand(2,1);
        
        O = molecules(1);
        N = molecules(2);
        CO = molecules(3);
        CN = molecules(4);
        C3N = molecules(5);
        NC2N = molecules(6);
        NC4N = molecules(7);
        OCCO = molecules(8);
        NCCO2 = molecules(9);
        
        %'CN numbers', CN
        %'k1', k_values(1)
        prop(1) = CN * (CN - 1) * k_values(1);
       % 'prop1', prop(1)
        
       %%%%%%%%% ESTABLISH CATALYSIS CONDITIONS FOR REACTION 2 %%%%%%%%%%%
       
        k_values(2) = uncatalysed;
        
        if C3N > 0 && NC2N > 0
            k_values(2) = catalysed * C3N * NC2N;
            
        elseif C3N > 0
            k_values(2) = catalysed * C3N;
        
        elseif NC2N > 0
            k_values(2) = catalysed * NC2N;
        end
        
        prop(2) = CN * C3N * k_values(2); %This reaction is catalysed by C3N and NC2N. kinetic constant is multiplied by catalyst numbers, as in Hordijk et al. (2014).
        prop(2);
        prop(3) = CN * O * CO * N * k_values(3);
        
       %%%%%%%%% ESTABLISH CATALYSIS CONDITIONS FOR REACTION 4 %%%%%%%%%%%

        
        k_values(4) = uncatalysed;
        
        if NC4N > 0
            k_values(4) = catalysed * NC4N;
        end
        
        prop(4) = CO * (CO - 1) * k_values(4); % Catalysed by NC4N.
        prop(5) = k_values(5);
        prop(6) = k_values(6);
        prop(7) = k_values(7);
        prop(8) = k_values(8);
        prop(9) = k_values(9);
        
        prop(10) = NC4N * k_values(10); %Lets enable these molecules to
        %accumulate.... so comment out the removal reactions!
        prop(11) = NC2N * k_values(11);
        prop(12) = N * k_values(12);
        
        prop(13) = NCCO2 * k_values(13);
        prop(14) = OCCO * k_values(14);
        
        
        alpha0 = sum(prop);
        %'Sum of propensities:', alpha0
        tau = (1/alpha0)*log(1/rr(1));
        time = time + tau;
        
        CumulativePropSum = cumsum(prop);
        
        if rr(2) < CumulativePropSum(1)
            %'rr(2) LESS THAN cumsum !'
            %'Executing reaction', 1
            molecules = molecules + (stoichiometric_matrix(:,1))'; %Choose first reaction if rr(2) < prop(1)/alpha0 !
        
        else
            for j = 1:length(prop)
                if rr(2) * alpha0 >= CumulativePropSum(j) && j < length(prop)
                    %'rr(2) GREATER OR EQUAL TO cumsum !'
                    if rr(2) * alpha0 < CumulativePropSum(j+1)
                        %'Executing reaction', j
                        %'Before the next reaction:', molecules
                     % disp(['Reaction' num2str(j+1) 'has occurred:'])
                        molecules = molecules + (stoichiometric_matrix(:,j+1))';
                    %molecules
                    %(stoichiometric_matrix(:,j+1))';
                    %'Propensities:', prop
                    end
                
                end
                
                
            end
        end
        
        k = k + 1;
        
        Oplot(k,i) = O;
        Nplot(k,i) = N;
        COplot(k,i) = CO;
        CNplot(k,i) = CN;
        C3Nplot(k,i) = C3N;
        NC2Nplot(k,i)= NC2N;
        NC4Nplot(k,i) = NC4N;
        OCCOplot(k,i) = OCCO;
        NCCO2plot(k,i) = NCCO2;
        
        timeplot(k,i) = time;
    end
                    
                    
end

figure(1);
set(gca,'Fontsize',20);
%h=stairs(timeplot(:,1),Oplot(:,1));

%set(h,'Color','r','Linewidth',1);
%h=stairs(timeplot(:,1),Nplot(:,1));
%set(h,'Color','m','Linewidth',1);
%h=stairs(timeplot(:,1),COplot(:,1));
%set(h,'Color','b','Linewidth',1);
%h=stairs(timeplot(:,1),CNplot(:,1));
%set(h,'Color','g','Linewidth',1);
%h=stairs(timeplot(:,1),C3Nplot(:,1));
%set(h,'Color','c','Linewidth',1);
h=stairs(timeplot(:,1),NC2Nplot(:,1));
set(h,'Color','k','Linewidth',1);
hold on
h=stairs(timeplot(:,1),NC4Nplot(:,1));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,1),OCCOplot(:,1));
set(h,'Color',[1,0.4,0.6],'Linewidth',1);
h=stairs(timeplot(:,1),NCCO2plot(:,1));
set(h,'Color',[1,0.8,0.6],'Linewidth',1);

'Number of non zero elements in NC2N vector:', nnz(NC2Nplot(:,1))

%NC4Nplot
xlabel('time [sec]');
ylabel('Number of Molecules');
legend('NC2N', 'NC4N', 'OCCO', 'NCCO2') % 'C3N', 'NC2N', 'NC4N', 'OCCO', 'NCCO2');
%axis([0 100 0 25]);
