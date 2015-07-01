function AutomatedGillespie

close all;

stoich_matrix = dlmread('stoichiometric_matrix.txt');

food_Numbers = dlmread('Food.txt');
non_food_Numbers = dlmread('non_food_molecules.txt');
catalysis_matrix = dlmread('catalysis_matrix.txt')
%global catalysis_matrix

for i = 1:length(food_Numbers)
    an_input_reaction = zeros(size(stoich_matrix, 1), 1);
    an_input_reaction(food_Numbers(i) + 1) = 1;
    
    stoich_matrix = [stoich_matrix an_input_reaction];
end

for i = 1:length(non_food_Numbers)
    an_output_reaction = zeros(size(stoich_matrix, 1), 1);
    an_output_reaction(non_food_Numbers(i) + 1) = -1;
    
    stoich_matrix = [stoich_matrix an_output_reaction];
end

stoich_matrix
Initial_k_values = zeros(1, size(stoich_matrix, 2));

uncatalysed = 0.1; 
catalysed = 1.0; %standard:= 1.0 s^-1
in = 1.2; % standard:= 0.3 s^-1
out = 0.01;
starting_number = 20;
loops = 50000;

for i = 1:size(stoich_matrix, 2) %for every reaction...
    NZelements = find(stoich_matrix(:,i));
    if length(NZelements) > 1
        Initial_k_values(i) = uncatalysed;
    elseif length(NZelements) == 1
        if stoich_matrix(NZelements, i) == 1
            Initial_k_values(i) = in;
        elseif stoich_matrix(NZelements, i) == -1
            Initial_k_values(i) = out;
        end
    else
        Initial_k_values(i) = 0; %If the reaction has no net product or reactant numbers!
    end
end

Initial_k_values

Initial_molNumbers = zeros(1, size(stoich_matrix, 1));

for i = 1:length(Initial_molNumbers)
    if ismember(i, food_Numbers + 1)
        Initial_molNumbers(i) = starting_number;
    end
end

Initial_molNumbers;
    
realisations = 5;
dist_colors = distinguishable_colors(length(non_food_Numbers));

for w = 1:realisations
 
    molNumbers = Initial_molNumbers;
    k_values = Initial_k_values; 
    time = 0;
    k = 1;
    
    for j = 1:length(molNumbers)
        Allplots(1, j, w) = molNumbers(j);
    end
    
    'allplots', Allplots(1,:,w);
    
    timeplot(1,w) = 0; %Now establish catalysis conditions, and update propensity functions.
    
    product_of_reactant_molecules = zeros(1, size(stoich_matrix, 2));
    propensities = zeros(1, size(stoich_matrix, 2));
     
    while k < loops
         
        rr = rand(2,1);
         
        for i = 1:size(catalysis_matrix, 2)
            catalysts = find(catalysis_matrix(:,i)); %Find catalyst indices
            if length(catalysts) >= 1
                
                for j = 1:length(catalysts)
                    catalyst_molecules = molNumbers(catalysts(j));
                    if catalyst_molecules > 0 
                        k_values(i) = catalysed * catalyst_molecules;
                    end
                end
            end
        end
        
        for i = 1:size(stoich_matrix, 2)
            reactants = find(stoich_matrix(:,i) < 0); %find reactant indices
            if length(reactants) > 0
                %'reactants', reactants;
                product_of_reactant_molecules(i) = 1;
                for j = 1:size(reactants)
                    reactant_stoichiometry = abs(stoich_matrix(reactants(j), i));
                    
                    if molNumbers(reactants(j)) >= reactant_stoichiometry
                        %factorial(molNumbers(reactants(j)) - abs(stoich_matrix(reactants(j), i)))
                        product_of_reactant_molecules(i) = product_of_reactant_molecules(i) * factorial(molNumbers(reactants(j)))/(factorial(molNumbers(reactants(j)) - reactant_stoichiometry));
                    else
                        product_of_reactant_molecules(i) = 0;
                    end
                end
            end
        end
        
        for i = 1:size(stoich_matrix, 2)
            if k_values(i) ~= in %Don't update the propensities for input reactions.
                propensities(i) = k_values(i) * product_of_reactant_molecules(i);
            else
                propensities(i) = k_values(i);
            end
        end
        
        alpha0 = sum(propensities);
        %'Sum of propensities:', alpha0
        tau = (1/alpha0)*log(1/rr(1));
        time = time + tau;
        
        CumulativePropSum = cumsum(propensities);
        
        if rr(2) < CumulativePropSum(1)
            molNumbers = molNumbers + (stoich_matrix(:,1))'; %Choose first reaction if rr(2) < prop(1)/alpha0 !
        
        else
            for j = 1:length(propensities)
                if rr(2) * alpha0 >= CumulativePropSum(j) && j < length(propensities)
                    if rr(2) * alpha0 < CumulativePropSum(j + 1)
                        %'Executing reaction', j
                        %'Before the next reaction:', molecules
                     % disp(['Reaction' num2str(j+1) 'has occurred:'])
                        molNumbers = molNumbers + (stoich_matrix(:,j + 1))';
                    %molecules
                    %(stoichiometric_matrix(:,j+1))';
                    %'Propensities:', prop
                    end
                
                end
                
                
            end
        end
        
        k = k + 1;
        
        for j = 1:length(molNumbers)
            Allplots(k, j, w) = molNumbers(j); %k-th loop, j-th molecule, w-th realisation.
        end
        
        timeplot(k, w) = time;
        
    end

%k_values
%product_of_reactant_molecules
%propensities
%molNumbers
    
    figure(w);
    set(gca,'Fontsize', 20);
    %a = 1;
    wanted_entries = [10 12 31 32];
    for j = 1:length(non_food_Numbers)
        'mol number:', non_food_Numbers(j)
        
        %if ismember(non_food_Numbers(j), wanted_entries)
        %    wanted_handles(a) = stairs(timeplot(:, w),Allplots(:, non_food_Numbers(j) + 1, w));
        %    set(wanted_handles(a),'Color',dist_colors(j,:),'Linewidth',1);
        %    a = a + 1;
        %    hold on
        %else
        %if ismember(non_food_Numbers(j) + 1, wanted_entries)
        h=stairs(timeplot(:, w),Allplots(:, non_food_Numbers(j) + 1, w)); %w = 1 before.
        set(h,'Color',dist_colors(j,:),'Linewidth',1);
        %if ismember(non_food_Numbers(j) + 1, wanted_entries)
        legendInfo{j} = ['M' num2str(non_food_Numbers(j))]; % or whatever is appropriate
            %a = a + 1;
        hold on
        
        %end
    end
    
    xlabel('Time [s]');
    ylabel('Number of Molecules');
    %legend(wanted_handles, 'M9', 'M11', 'M30', 'M31');  
    legend(legendInfo)
    
end

%%%DETERMINISTIC SOLUTION FOR ITER=1 TOY-CHEMICAL REACTION SYSTEM%%%

% [tdet,det] = ode45(@myode, [0 500], Initial_molNumbers, [], Initial_k_values, stoich_matrix);
% 
% figure(6);
% set(gca,'Fontsize',20);
% 
% plot(tdet,det(:,2),'--k','Linewidth',4); %plot molecule 2!
% hold on;
% 
% h=stairs(timeplot(:,1),Allplots(:, 2, 1));
% set(h,'Color','r','Linewidth',1);
% h=stairs(timeplot(:,2),Allplots(:, 2, 2));
% set(h,'Color','m','Linewidth',1);
% h=stairs(timeplot(:,3),Allplots(:, 2, 3));
% set(h,'Color','b','Linewidth',1);
% h=stairs(timeplot(:,4),Allplots(:, 2, 4));
% set(h,'Color','g','Linewidth',1);
% h=stairs(timeplot(:,5),Allplots(:, 2, 5));
% set(h,'Color','c','Linewidth',1);
% 
% 
% plot(tdet,det(:,2),'--k','Linewidth',4);
% xlabel('Time [sec]');
% ylabel('Number of OCC=O Molecules');
% legend('solution of ODEs');
% 
% figure(7);
% set(gca,'Fontsize',20);
% plot(tdet,det(:,5),'--k','Linewidth',4); %plot molecule 5!
% hold on;
% 
% h=stairs(timeplot(:,1),Allplots(:, 5, 1));
% set(h,'Color','r','Linewidth',1);
% h=stairs(timeplot(:,2),Allplots(:, 5, 2));
% set(h,'Color','m','Linewidth',1);
% h=stairs(timeplot(:,3),Allplots(:, 5, 3));
% set(h,'Color','b','Linewidth',1);
% h=stairs(timeplot(:,4),Allplots(:, 5, 4));
% set(h,'Color','g','Linewidth',1);
% h=stairs(timeplot(:,5),Allplots(:, 5, 5));
% set(h,'Color','c','Linewidth',1);
% 
% plot(tdet,det(:,5),'--k','Linewidth',4);
% xlabel('Time [sec]');
% ylabel('Number of N=C(C#C)C#N Molecules');
% legend('solution of ODEs');
% 
% 
% function dydt = myode(t, M, Initial_k_values, stoich_matrix)
% 
% %%%pass Initial_k_values and stoich_matrix to this function!!!%%%
% 
% k = Initial_k_values;
% 
% if M(7) > 0
%     k(2) = 1.0 * M(7);
% else
%     k(2) = 0.1;
% end
% 
% if M(5) > 0
%     k(1) = 1.0 * M(5);
% else
%     k(1) = 0.1;
% end
%        
%     
% 
% flux_vector = [k(1)*M(1)^2; k(2)*M(3)*M(4); k(3)*M(4)^2; k(4)*M(1)*M(4)*M(7); k(5); k(6); k(7); k(8); k(9); k(10)*M(5); k(11)*M(6); k(12)*M(8); k(13)*M(2)];
% 
% dydt = stoich_matrix * flux_vector;


%%%DETERMINISTIC SOLUTION FOR ITER=2 TOY-CHEMICAL REACTION SYSTEM%%%

[tdet,det] = ode15s(@myode, [0 500], Initial_molNumbers, [], Initial_k_values, stoich_matrix);

figure(6);
set(gca,'Fontsize',20);

plot(tdet,det(:,15),'--k','Linewidth',4); %plot molecule 14!
hold on;

h=stairs(timeplot(:,1),Allplots(:, 15, 1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Allplots(:, 15, 2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Allplots(:, 15, 3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Allplots(:, 15, 4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Allplots(:, 15, 5));
set(h,'Color','c','Linewidth',1);


%plot(tdet,det(:,2),'--k','Linewidth',4);
xlabel('Time [sec]');
ylabel('Number of NCC(O)=O Molecules');
legend('solution of ODEs');

figure(7);
set(gca,'Fontsize',20);
plot(tdet,det(:,10),'--k','Linewidth',4); %plot molecule 9!
hold on;

h=stairs(timeplot(:,1),Allplots(:, 10, 1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Allplots(:, 10, 2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Allplots(:, 10, 3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Allplots(:, 10, 4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Allplots(:, 10, 5));
set(h,'Color','c','Linewidth',1);

%plot(tdet,det(:,5),'--k','Linewidth',4);
xlabel('Time [sec]');
ylabel('Number of NCC#N Molecules');
legend('solution of ODEs');

figure(8);
set(gca,'Fontsize',20);
plot(tdet,det(:,3),'--k','Linewidth',4); %plot molecule 2!
hold on;

h=stairs(timeplot(:,1),Allplots(:, 3, 1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Allplots(:, 3, 2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Allplots(:, 3, 3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Allplots(:, 3, 4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Allplots(:, 3, 5));
set(h,'Color','c','Linewidth',1);

%plot(tdet,det(:,5),'--k','Linewidth',4);
xlabel('Time [sec]');
ylabel('Number of OCC(O)C=O Molecules');
legend('solution of ODEs');

figure(9);
set(gca,'Fontsize',20);
plot(tdet,det(:,11),'--k','Linewidth',4); %plot molecule 10!
hold on;

h=stairs(timeplot(:,1),Allplots(:, 11, 1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Allplots(:, 11, 2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Allplots(:, 11, 3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Allplots(:, 11, 4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Allplots(:, 11, 5));
set(h,'Color','c','Linewidth',1);

%plot(tdet,det(:,5),'--k','Linewidth',4);
xlabel('Time [sec]');
ylabel('Number of NCC(=N)C#N Molecules');
legend('solution of ODEs');

figure(10);
set(gca,'Fontsize',20);
plot(tdet,det(:,7),'--k','Linewidth',4); %plot molecule 6!
hold on;

h=stairs(timeplot(:,1),Allplots(:, 7, 1));
set(h,'Color','r','Linewidth',1);
h=stairs(timeplot(:,2),Allplots(:, 7, 2));
set(h,'Color','m','Linewidth',1);
h=stairs(timeplot(:,3),Allplots(:, 7, 3));
set(h,'Color','b','Linewidth',1);
h=stairs(timeplot(:,4),Allplots(:, 7, 4));
set(h,'Color','g','Linewidth',1);
h=stairs(timeplot(:,5),Allplots(:, 7, 5));
set(h,'Color','c','Linewidth',1);

%plot(tdet,det(:,5),'--k','Linewidth',4);
xlabel('Time [sec]');
ylabel('Number of N=CC#N Molecules');
legend('solution of ODEs');


function dydt = myode(t, M, Initial_k_values, stoich_matrix)

k = Initial_k_values;

%%%CATALYSIS CONDITIONS FOR THE PARTICULAR CRS%%%

if M(10) > 0 && M(11) > 0 && M(13) > 0
    k(8) = 1.0 * M(10) * M(11) * M(13);
    
elseif M(10) > 0 && M(11) > 0
    k(8) = 1.0 * M(10) * M(11);
    
elseif M(11) > 0 && M(13) > 0
    k(8) = 1.0 * M(11) * M(13);
    
elseif M(10) > 0 && M(13) > 0
    k(8) = 1.0 * M(10) * M(13);
elseif M(10) > 0
    k(8) = 1.0 * M(10);
elseif M(11) > 0
    k(8) = 1.0 * M(11);
elseif M(13) > 0
    k(8) = 1.0 * M(13);
end
    
if M(7) > 0 && M(11) > 0
    k(9) = 1.0 * M(7) * M(11);
    k(1) = 1.0 * M(7);
elseif M(7) > 0
    k(9) = 1.0 * M(7);
    k(1) = 1.0 * M(7);
elseif M(11) > 0
    k(9) = 1.0 * M(11);
end

if M(4) > 0 && M(5) > 0
    k(6) = 1.0 * M(4) * M(5);
    k(5) = 1.0 * M(4);
elseif M(4) > 0
    k(6) = 1.0 * M(4);
    k(5) = 1.0 * M(4);
elseif M(5) > 0
    k(6) = 1.0 * M(5);
end

if M(13) > 0
    k(11) = 1.0 * M(13);
end
%global catalysis_matrix;

%k = M * catalysis_matrix;

flux_vector = [k(1)*(M(1))^2; k(2)*M(1)*M(2); k(3); k(4); k(5)*M(4)*M(5); k(6)*(M(5))^2; k(7)*M(5)*M(6); k(8)*M(5)*M(7); k(9)*M(5)*M(10); k(10)*M(1)*M(5)*M(12); k(11)*M(2)*M(5)*M(12); k(12)*M(10)*(M(13))^2; k(13); k(14); k(15); k(16); k(17); k(18)*M(6); k(19)*M(8); k(20)*M(7); k(21)*M(9); k(22)*M(14); k(23)*M(10); k(24)*M(11); k(25)*M(15); k(26)*M(3); k(27)*M(2)];

dydt = stoich_matrix * flux_vector;

%%%DETERMINISTIC SOLUTION FOR ITER=3 TOY-CHEMICAL REACTION SYSTEM%%%

% [tdet,det] = ode45(@myode, [0 500], Initial_molNumbers, [], Initial_k_values, stoich_matrix);
% 
% figure(6);
% set(gca,'Fontsize',20);
% 
% plot(tdet,det(:,8),'--k','Linewidth',4); %plot molecule 7!
% hold on;
% 
% h=stairs(timeplot(:,1),Allplots(:, 8, 1));
% set(h,'Color','r','Linewidth',1);
% h=stairs(timeplot(:,2),Allplots(:, 8, 2));
% set(h,'Color','m','Linewidth',1);
% h=stairs(timeplot(:,3),Allplots(:, 8, 3));
% set(h,'Color','b','Linewidth',1);
% h=stairs(timeplot(:,4),Allplots(:, 8, 4));
% set(h,'Color','g','Linewidth',1);
% h=stairs(timeplot(:,5),Allplots(:, 8, 5));
% set(h,'Color','c','Linewidth',1);
% 
% 
% %plot(tdet,det(:,2),'--k','Linewidth',4);
% xlabel('Time [sec]');
% ylabel('Number of N=C(C#C)C#N Molecules');
% legend('solution of ODEs');
% 
% figure(7);
% set(gca,'Fontsize',20);
% plot(tdet,det(:,23),'--k','Linewidth',4); %plot molecule 22!
% hold on;
% 
% h=stairs(timeplot(:,1),Allplots(:, 23, 1));
% set(h,'Color','r','Linewidth',1);
% h=stairs(timeplot(:,2),Allplots(:, 23, 2));
% set(h,'Color','m','Linewidth',1);
% h=stairs(timeplot(:,3),Allplots(:, 23, 3));
% set(h,'Color','b','Linewidth',1);
% h=stairs(timeplot(:,4),Allplots(:, 23, 4));
% set(h,'Color','g','Linewidth',1);
% h=stairs(timeplot(:,5),Allplots(:, 23, 5));
% set(h,'Color','c','Linewidth',1);
% 
% %plot(tdet,det(:,5),'--k','Linewidth',4);
% xlabel('Time [sec]');
% ylabel('Number of NCC(O)=O Molecules');
% legend('solution of ODEs');
% 
% figure(8);
% set(gca,'Fontsize',20);
% plot(tdet,det(:,22),'--k','Linewidth',4); %plot molecule 21!
% hold on;
% 
% h=stairs(timeplot(:,1),Allplots(:, 22, 1));
% set(h,'Color','r','Linewidth',1);
% h=stairs(timeplot(:,2),Allplots(:, 22, 2));
% set(h,'Color','m','Linewidth',1);
% h=stairs(timeplot(:,3),Allplots(:, 22, 3));
% set(h,'Color','b','Linewidth',1);
% h=stairs(timeplot(:,4),Allplots(:, 22, 4));
% set(h,'Color','g','Linewidth',1);
% h=stairs(timeplot(:,5),Allplots(:, 22, 5));
% set(h,'Color','c','Linewidth',1);
% 
% %plot(tdet,det(:,5),'--k','Linewidth',4);
% xlabel('Time [sec]');
% ylabel('Number of NC(CO)C(O)=O Molecules');
% legend('solution of ODEs');
% 
% figure(9);
% set(gca,'Fontsize',20);
% plot(tdet,det(:,5),'--k','Linewidth',4); %plot molecule 4!
% hold on;
% 
% h=stairs(timeplot(:,1),Allplots(:, 5, 1));
% set(h,'Color','r','Linewidth',1);
% h=stairs(timeplot(:,2),Allplots(:, 5, 2));
% set(h,'Color','m','Linewidth',1);
% h=stairs(timeplot(:,3),Allplots(:, 5, 3));
% set(h,'Color','b','Linewidth',1);
% h=stairs(timeplot(:,4),Allplots(:, 5, 4));
% set(h,'Color','g','Linewidth',1);
% h=stairs(timeplot(:,5),Allplots(:, 5, 5));
% set(h,'Color','c','Linewidth',1);
% 
% %plot(tdet,det(:,5),'--k','Linewidth',4);
% xlabel('Time [sec]');
% ylabel('Number of OCC(=O)CO Molecules');
% legend('solution of ODEs');
% 
% function dydt = myode(t, M, Initial_k_values, stoich_matrix)
% 
% k = Initial_k_values;
% 
% %%%CATALYSIS CONDITIONS FOR THE PARTICULAR CRS%%%
% 
% if M(1) > 0 && M(10) > 0
%     k(2) = 1.0 * M(1) * M(10);
% elseif M(1) > 0
%     k(2) = 1.0 * M(1);
% elseif M(10) > 0
%     k(2) = 1.0 * M(10);
% end
% 
% if M(1) > 0 && M(18) > 0
%     k(4) = 1.0 * M(1) * M(18);
% elseif M(1) > 0
%     k(4) = 1.0 * M(1);
% elseif M(18) > 0
%     k(4) = 1.0 * M(18);
% end
% 
% if M(9) > 0 && M(12) > 0 && M(21) > 0
%     k(3) = 1.0 * M(9) * M(12) * M(21);
% elseif M(9) > 0 && M(12) > 0
%     k(3) = 1.0 * M(9) * M(12);
% elseif M(9) > 0 && M(21) > 0
%     k(3) = 1.0 * M(9) * M(21);
% elseif M(21) > 0 && M(12) > 0
%     k(3) = 1.0 * M(21) * M(12);
% elseif M(9) > 0
%     k(3) = 1.0 * M(9);
% elseif M(12) > 0
%     k(3) = 1.0 * M(12);
% elseif M(21) > 0
%     k(3) = 1.0 * M(21);
% end
% 
% if M(18) > 0
%     k(6) = 1.0 * M(18);
% end
% 
% if M(3) > 0
%     k(8) = 1.0 * M(3);
% end
% 
% if M(8) > 0
%     k(10) = 1.0 *M(8);
% end
% 
% if M(20) > 0 && M(21) > 0
%     k(11) = 1.0 * M(20) * M(21);
% elseif M(20) > 0
%     k(11) = 1.0 * M(20);
% elseif M(21) > 0
%     k(11) = 1.0 * M(21);
% end
% 
% if M(9) > 0 && M(13) > 0
%     k(12) = 1.0 * M(9) * M(13);
% elseif M(9) > 0
%     k(12) = 1.0 * M(9);
% elseif M(13) > 0
%     k(12) = 1.0 * M(13);
% end
% 
% if M(2) > 0 && M(7) > 0
%     k(13) = 1.0 * M(2) * M(7);
% elseif M(2) > 0
%     k(13) = 1.0 * M(2);
% elseif M(7) > 0
%     k(13) = 1.0 * M(7);
% end
% 
% if M(18) > 0
%     k(15) = 1.0 * M(18);
% end
% 
% catalysts_reaction16 = [M(1) M(4) M(5) M(10)];
% catalysts_present = find(catalysts_reaction16 > 0);
% if length(catalysts_present) > 0
%     k(16) = 1.0;
%     for i = 1:length(catalysts_present)
%         k(16) = k(16) * catalysts_reaction16(catalysts_present(i));
%     end
% end
% 
% if M(19) > 0
%     k(17) = 1.0 * M(19);
% end
% 
% if M(16) > 0
%     k(18) = 1.0 * M(16);
% end
% 
% if M(4) > 0 && M(14) > 0
%     k(19) = 1.0 * M(4) * M(14);
% elseif M(4) > 0
%     k(19) = 1.0 * M(4);
% elseif M(14) > 0
%     k(19) = 1.0 * M(14);
% end
% 
% 
% 
% flux_vector = [k(1)*M(1)^2; k(2)*M(1)*M(3); k(3)*M(1)*M(2); 0; k(5)*M(3); 0; k(7)*M(6)*M(7); k(8)*M(7)^2; k(9)*M(7)*M(8); k(10)*M(7)*M(10); k(11)*M(7)*M(9); k(12)*M(7)*M(12); k(13)*M(7)*M(14); k(14)*M(7)*M(16); k(15)*M(7)*M(17); k(16)*M(1)*M(7)*M(19); k(17)*M(3)*M(7)*M(19); k(18)*M(2)*M(7)*M(19); k(19)*M(14)*M(20)^2; k(20)*M(16)*M(20)^2; k(21); k(22); k(23); k(24); k(25); k(26)*M(8); k(27)*M(10); k(28)*M(11); k(29)*M(9); k(30)*M(12); k(31)*M(13); k(32)*M(21); k(33)*M(14); k(34)*M(15); k(35)*M(22); k(36)*M(16); k(37)*M(17); k(38)*M(18); k(39)*M(23); k(40)*M(5); k(41)*M(4); k(42)*M(3); k(43)*M(2)];
% 
% dydt = stoich_matrix * flux_vector;

