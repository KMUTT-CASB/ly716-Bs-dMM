function low_concentration_extracellular_met_model
global tdat Ydat

% Function to model low concentration

% Initialize initial value from experiment
Y4(1) = 0.1;
Y4(2) = 292.397660818713;
Y4(3) = 0;
Y4(4) = 0;
Y4(5) = 0;
Y4(6) = 0;

% Parameters from 100 g/L

% Using Monod
k(1) =  0.75;
k(2) =  0.7;
k(3) =  7.3;
k(4) = 18.19;
k(5) =  6.417;
k(6) = 1546.455;
k(7) =  0.252;
k(8) =  63;
k(11) = 0.031;
k(14) = 0.25;
k(15) =  1533.588;
k(16) = 0.096/2;

% Function to model dynamic extracellular metabolites

    function dC4=DifEq4(t,Y4)
            dc4dt=zeros(6,1);
            u= (0.2454.*Y4(2))./(Y4(2)+ k(1)); % Monod Equation
            p2 = 0.096*Y4(2); % Hydrolysis reaction, levansucrase
            p3= ((k(5).*(Y4(2))./((Y4(2)) + (k(6)))).*k(4)); % Transfructosylation reaction, levansucrase
            p4= (k(7).*((Y4(5)-14.798))./(((Y4(5)-14.798)) + (k(8)))).*k(4); % Levan degradation, levansucrase
           
    dc4dt(1) = (u-0.01).*Y4(1); % Biomass Growth with death constant
    dc4dt(2) = -(((u)*Y4(1))/0.025)-0.5*p2-(((p3)*Y4(1))); % Sucrose cons
    dc4dt(3) = (p2*rectangularPulse(0,6,t))+(((p3)*Y4(1))); % Glucose prod
    dc4dt(4) = (p2*rectangularPulse(0,6,t))+((p4)*Y4(1)*heaviside(t-18)); % Fructose prod
    dc4dt(5) = ((p3)*Y4(1))-(p4*Y4(1)*heaviside(t-18)); % Levan Production
    dc4dt(6) = u*k(4)*Y4(1); % Levansucrase over time based on biomass growth
    dC4=dc4dt;
    end

% Solving equation
tdat = 0:1:120;
[t,Y4] = ode23(@DifEq4,tdat,[Y4(1) Y4(2) Y4(3) Y4(4) Y4(5) Y4(6)]);

% Import dataset
data = xlsread('low_concentration_extracellular_met.xlsx',1);
tdat = data(:,1);
Ydat = [data(:,2) data(:,3) data(:,4) data(:,5) data(:,6) data(:,7) data(:,8)];

% Set gray color
grayColor = '#808080';

figure(2)
plot(t(:,1),Y4(:,1),'Color',grayColor,'LineWidth',2,'MarkerSize',8)
hold on
errorbar(tdat(:,1),Ydat(:,1),Ydat(:,6),'Color',grayColor,'Marker','o','LineStyle','none','LineWidth',2,'MarkerSize',8)
ylim([0 3])
xlim([0 120])
xticks(0:30:120)
yticks(0:1:3)
xlabel('Time (h)')
ylabel('Biomass (gDW L^{-1})')
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(3)
plot(t(:,1),Y4(:,2),'Color',grayColor,'LineWidth',2,'MarkerSize',8)
hold on
plot(tdat(:,1),Ydat(:,2),'o','Color',grayColor,'LineWidth',2,'MarkerSize',8)
xlabel('Time (h)')
ylabel('Sucrose (mmol L^{-1})')
ylim([0 300])
xlim([0 120])
xticks(0:30:120)
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(4)
plot(t(:,1),Y4(:,3),'Color',grayColor,'LineWidth',2,'MarkerSize',8)
hold on
plot(tdat(:,1),Ydat(:,3),'o','Color',grayColor,'LineWidth',2,'MarkerSize',8)
xlabel('Time (h)')
ylabel('Glucose (mmol L^{-1})')
ylim([0 250])
xlim([0 120])
xticks(0:30:120)
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(5)
plot(t(:,1),Y4(:,5),'Color',grayColor,'LineWidth',2,'MarkerSize',8)
hold on
errorbar(tdat(:,1),Ydat(:,5),Ydat(:,7),'Color',grayColor,'Marker','o','LineStyle','none','LineWidth',2,'MarkerSize',8)
xlabel('Time (h)')
ylabel('Levan (mmol L^{-1})')
ylim([0 80])
xlim([0 120])
xticks(0:30:120)
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(6)
plot(t(:,1),Y4(:,4),'Color',grayColor,'LineWidth',2,'MarkerSize',8)
hold on
plot(tdat(:,1),Ydat(:,4),'o','Color',grayColor,...
    'LineWidth',2,'MarkerSize',8)
xlabel('Time (h)')
ylabel('Fructose (mmol L^{-1})')
ylim([0 250])
xlim([0 120])
xticks(0:30:120)
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

figure(7)
plot(t(:,1),Y4(:,6),'Color',grayColor,'LineWidth',2,'MarkerSize',8,...
    'LineWidth',2)
xlabel('Time (h)')
ylabel({'Simulated Enzyme';'(mg levansucrase L^{-1})'})
ylim([0 60])
xlim([0 120])
xticks(0:30:120)
set(gca, 'box', 'off')
fontname("Arial")
fontsize(24,"points")

end