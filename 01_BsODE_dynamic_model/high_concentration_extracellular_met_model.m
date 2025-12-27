function high_concentration_extracellular_met_model
global tdat Ydat

% Function to model high concentration

% Initialize initial value from experiment
Y5(1) = 0.1;
Y5(2) = 730.994152;
Y5(3) = 0;
Y5(4) = 0;
Y5(5) = 0;
Y5(6) = 0;

% Parameters from 250 g/L

% Using Monod
k(1) =  0.75;
k(2) =  0.096; 
k(3) = 45;
k(4) = 6.417;
k(5) =  1546.45;
k(6) = 0.167;
k(7) =  1533.588;
k(8) =  0.252;
k(9) =  63;
k(10) = 0.006;

% Function to model dynamic extracellular metabolites

function dC5=DifEq5(t,Y5)
            dc5dt=zeros(6,1);
            u= (0.11*Y5(2))./(Y5(2)+ k(1)); % Monod Equation
            p2 = k(2)*Y5(2); % Hydrolysis reaction
            p3= (((k(4)*(Y5(2))./k(5))+(((k(6)*((Y5(4))*(Y5(2))))./(k(7)))))...
                ./(1 + ((Y5(2))./k(5)) + ((Y5(4))*(Y5(2)))./(k(7))))*k(3); % Multiple binding sites for levan synthesis, sucrose has higher affinity
            p4= (k(8)*((Y5(5)-67.394))./(((Y5(5)-67.394)) + (k(9))))*(47.042/4); % Levan degradation
           

    dc5dt(1) = (u-0.01)*Y5(1); % Biomass Growth with death constant
    dc5dt(2) = -((((u)*Y5(1))/k(10)))-(((p3)*Y5(1))); % Sucrose cons
    dc5dt(3) = (((p2*rectangularPulse(0,6,t))))+((p3)*Y5(1)); % Glucose pro
    dc5dt(4) = (((p2*rectangularPulse(0,6,t))))-((p3)*Y5(1))+((p4)*heaviside(t-36)); %Fructose prod
    dc5dt(5) = ((p3)*Y5(1))-(p4*heaviside(t-36)); % Levan Production
    dc5dt(6) = u*k(3)*Y5(1); % Levansucrase over time based on biomass growth
    dC5=dc5dt;
        end

% Solving equation
tdat = 0:1:120;
[t,Y5] = ode45(@DifEq5,tdat,[Y5(1) Y5(2) Y5(3) Y5(4) Y5(5) Y5(6)]);

% Import dataset
data = xlsread('high_concentration_extracellular_met.xlsx',1);
tdat = data(:,1);
Ydat = [data(:,2) data(:,3) data(:,4) data(:,5) data(:,6) data(:,7) data(:,8)];

% Plot the dynamic model

figure(2)
plot(t(:,1),Y5(:,1),'black-','LineWidth',2,'MarkerSize',8)
hold on
errorbar(tdat(:,1),Ydat(:,1),Ydat(:,6),'blacko','LineStyle','none','LineWidth',2,'MarkerSize',8)
ylim([0 4])
xlim([0 120])
xticks(0:30:120)
yticks(0:1:4)
xlabel('Time (h)')
ylabel('Biomass (gDW L^{-1})')
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(3)
plot(t(:,1),Y5(:,2),'black-',...
         tdat(:,1),Ydat(:,2),'blacko','LineWidth',2,'MarkerSize',8)
xticks(0:30:120)
xlabel('Time (h)')
ylabel('Sucrose (mmol L^{-1})')
ylim([0 1000])
xlim([0 120])
legend('Simulation','Experiment')
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")


figure(4)
plot(t(:,1),Y5(:,3),'black-',...
      tdat(:,1),Ydat(:,3),'blacko','LineWidth',2,'MarkerSize',8)
xticks(0:30:120)
xlabel('Time (h)')
ylabel('Glucose (mmol L^{-1})')
ylim([0 1000])
xlim([0 120])
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(5)
plot(t(:,1),Y5(:,5),'black-','LineWidth',2,'MarkerSize',8)
hold on
errorbar(tdat(:,1),Ydat(:,5),Ydat(:,7),'blacko','LineStyle','none','LineWidth',2,'MarkerSize',8)
xticks(0:30:120)
xlabel('Time (h)')
ylabel('Levan (mmol L^{-1})')
ylim([0 300])
xlim([0 120])
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")


figure(6)
plot(t(:,1),Y5(:,4),'black-',...
     tdat(:,1),Ydat(:,4),'blacko','LineWidth',2,'MarkerSize',8)
ylim([0 500])
xticks(0:30:120)
xlabel('Time (h)')
ylabel('Fructose (mmol L^{-1})')
ylim([0 600])
xlim([0 120])
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(7)
plot(t(:,1),Y5(:,6),'black-',...
    'LineWidth',2)
ylim([0 200])
xlim([0 120])
xticks(0:30:120)
xlabel('Time (h)')
ylabel({'Simulated Enzyme';'(mg levansucrase L^{-1})'})
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")


end