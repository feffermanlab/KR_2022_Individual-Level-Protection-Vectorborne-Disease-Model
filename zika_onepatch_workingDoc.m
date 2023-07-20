% Testing expanded vector-control model - One patch
% -> includes personal protection (protected susceptibles w/ lower transmission)

clear 
close all

%%

% Set parameters (based on values from Suarez, 2020)
b = (9/1000)/365;
betaH = 1.5*10^-4;
muH = (8.6/1000)/365;
r = 0.037;
omega = 0;
delta = 5.468913e-05;
n = 10;
nu = 1/7;
betaM = 3*10^-4;
muM = 1/13;

% Variable parameters
rho = 0.2;
K = 20000 + randi(5000,1,1);
gamma = 0;
e = 0;
gamma_D = 1.5/700;
gamma_B = 0.1/1200;
lamdaP = 1/2;
lag = 0;    % tau

% Initial conditions
N = 700;
Sh0 = N - 1;
Sp0 = 0;
Ih0 = 1;
Ip0 = 0;
Rh0 = 0;
Dh0 = 0;
Lm0 = 0;
Sm0 = 1200;
Im0 = 0;
Cm0 = 0;
Cl0 = 0;
%Cp0 = 0;

IC = [Sh0 Sp0 Ih0 Ip0 Rh0 Dh0 Lm0 Sm0 Im0 Cm0 Cl0];

t_window = 1:365;

[t,x] = ode23s(@zika_model_nocontrol,t_window,IC,[],b,betaH,muH,rho,r,omega,delta,n,K,nu,betaM,muM,gamma,e,gamma_D,gamma_B,lamdaP,lag);

%% 

FigH = figure('Position', get(0, 'Screensize'));
hold on

subplot(1,2,1)
plot(t_window,x(:,1),'LineWidth',2)

hold on

plot(t_window,x(:,2),'LineWidth',2)
plot(t_window,x(:,3),'LineWidth',2)
plot(t_window,x(:,4),'LineWidth',2)

legend({'Fully susceptible','Protected susceptible','Infectious','Protected Infectious'})
title('Humans')
set(gca,'FontSize', 20);

xlim([0 length(t_window)])

subplot(1,2,2)
hold on

yyaxis left
line6 = plot(t_window,x(:,7),'--','LineWidth',1.5)
line7 = plot(t_window,x(:,8),'-','LineWidth',1.5)

yyaxis right
line8 = plot(t_window,x(:,9),'LineWidth',2)

legend([line6 line7 line8],{'Larvae','Susceptible adults','Infectious adults'})
title('Mosquitos')
set(gca,'FontSize', 20);

xlim([0 length(t_window)])
set(gcf,'Color','White')


%saveas(FigH, '')
