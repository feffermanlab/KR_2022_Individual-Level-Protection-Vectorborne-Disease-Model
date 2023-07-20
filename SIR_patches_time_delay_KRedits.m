function [S_out, I_out, R_out]  = SIR_patches_time_delay_KRedits(folder, information, MaxTime, p, alpha, gamma, epsilon, N, tau, seed, network)
% This function calculates and plots the temporal evolution of the model 
% The output is the values of human succeptible population and
% recovered population at the end of the simulation, 
% and  the peak of infected people for every patch, as a function of the
% environmental concern parameter (epsilon). 
% You can include a delay on the application of control measures changing
% the value of variable "tau"

%close all

% All the constants that I will not change
beta    = [0.15 0.30]/1000;       % transmision rate for [humans mosquitoes]
mu      = [2.35616e-05  1/13];    % natural death rate for [humans mosquitoes]
b       = 2.4657534e-05;          % human birth rate
r       = 0.037;                  % human recovery rate
w       = 0;                      % human death rate from disease
delta   = 5.468913e-05;           % composite rate 
nu      = 1/7;                    % maturation rate
eta     = 10;                     % egg laying rate - it says 10 on the table

m       = 12;                     % number of equations

if nargin == 0                          % if the number of input arguments is zero
    seed     = randi(intmax('int32'),1,1)
    MaxTime  = 500;                     % end of simulation
    %epsilon  = [150 150] ;   %High     % environment concern demotivates [mosquitos  larvae] pesticide
    %epsilon = [500 500];   % Low
    epsilon = [200 200];   %Medium
    gamma    = exp(-epsilon/50);        % motivation to control [mosquitos  larvae] due to Infected humans
    alpha    = 100 * gamma;             % motivation to control [mosquitos  larvae] due to Severe Outcomes (always 100 times larger than the regular infected)
    N        = 25;                      % number of patches
    tau      = 7;                       % time delay
    p        = 0.2;                     % fraction of people that moves to other patches
    
    %Additional parameters for expanded model
    rho      = 0.2;                     % Reduction in transmission for protected susceptibles 
    gamma_D  = 15/700;                 % Motivation to use personal protection based on fear of disease
    gamma_B  = 0.1/1200;                 % Motivation to use personal protection for fear of being bitten
    lamda    = 1/4;                     % 1/ length of use of personal protection
    
    network  = 'SQ';                    % type of network: SQ, SF, full
    information = 'global';             % 'global' or 'local' information
    folder = join([network, information, 'N', string(N),'p', string(p), 't', string(MaxTime)], {'-'} );
    mkdir(folder);
end

rng(seed) % seed for pseudorandom number generator

gB_short = trunc(gamma_B,6);
gD_short = trunc(gamma_D,6);

%filename = strcat(folder, '/2-a1-', string(alpha(1)), '-a2-', string(alpha(2)) ,'-g1-', string(gamma(1)),'-g2-', string(gamma(2)),'-e1-', string(epsilon(1)),'-e2-', string(epsilon(2)),'-delay-',string(tau),'-rho-',string(rho));
filename = strcat(folder,'MLL_lamdaEQ_4');

%initial conditions in each patch
X0 = [];
Kl = [];
%delay = [];
for i = 1:N
       %  [R0 appears two times, the first one is for recovered people as a
       %  function of time in a given patch, the second one is reserved to
       %  accumulate the total number of recovered people since the begining
       %  [S0 Humans           ,Sp0,Ih0,Ip0,R0, D0, L0, S0 mosq.,             I0, Cm0, Cl0, R0]
    ini = [700 + randi(100,1,1); 0; 0;  0;  0;  0;  0; 1200 + randi(500,1,1); 0;   0;   0;   0];
    X0  = [X0; ini];
    Kl  = [Kl; 20000 + randi(5000,1,1)];  %carrying capacity of each patch: 20000 + randi(5000,1,1)
end

delay = tau * ones(1,2*N);  % an array of ones times the delay. It is the same for all
                            % the 2*N variables of control efforts. Cm(j) and Cl(j) 

%%------------ Network and mobility --------------------------------------
% define the network structure and mobility conditions for different cases

if strcmp(network, 'SQ')
    lat_size = sqrt(N);
    mov = square_adj_mat(lat_size);
    %Infect one person in the middle of the lattice
    J0 = (lat_size^2 + 1) / 2;
    dd = distance(mov); % dd(i,j) = dd(j,i) is the distance from patch i to j
end

X0((J0-1)*m+3) = 1;  % J0 is the first infected patch


for jj=1:N
    mov(jj,:)=mov(jj,:)/sum(mov(jj,:));  %normalized by rows. (people leaving can't be more than 100%)
end 
mov = p * mov;    % p is the fraction of people that leave the patch
%%------------------------ END --------------------------------------------
                            
                            
                                %temporal span
%tspan = 0:MaxTime/500:MaxTime;  %500 points from 0 to MaxTime
tspan = 0:MaxTime;

%solving the equations
if strcmp(information, 'global')
    if tau == 0 
       options = odeset('RelTol',1e-3,'AbsTol',1e-3);
       pop = ode23s(@SIR_eq_no_delay_global, tspan, X0, options);
    else
        pop = dde23(@SIR_eq_delay_global, delay, X0, tspan);
    end
end


% separate the results in different variables
Sh=pop.y(1:m:m*N,:)';
Sp=pop.y(2:m:m*N,:)';                     
Ih=pop.y(3:m:m*N,:)';
Ip=pop.y(4:m:m*N,:)';
Rh=pop.y(12:m:m*N,:)';
Dh=pop.y(6:m:m*N,:)';
Lm=pop.y(7:m:m*N,:)';
Sm=pop.y(8:m:m*N,:)';
Im=pop.y(9:m:m*N,:)';
Cm=pop.y(10:m:m*N,:)';
Cl=pop.y(11:m:m*N,:)';


% % Normalize by the total number of people on each patch
N_total = zeros(length(Sh(:,1)),1);
M_total = zeros(length(Sh(:,1)),1);
for j = 1:N
    N_total(:) = Sh(:,j) + Sp(:,j) + Ih(:,j) + Ip(:,j) + Rh(:,j);
    M_total(:) = Sm(:,j) + Im(:,j);
    Sh(:,j) = Sh(:,j) ./ N_total;
    Sp(:,j) = Sp(:,j) ./ N_total;
    Ih(:,j) = Ih(:,j) ./ N_total;
    Ip(:,j) = Ip(:,j) ./ N_total;
    Rh(:,j) = Rh(:,j) ./ N_total;
    Dh(:,j) = Dh(:,j) ./ N_total;
    Lm(:,j) = Lm(:,j) ./ Kl(j);
    Sm(:,j) = Sm(:,j) ./ M_total;
    Im(:,j) = Im(:,j) ./ M_total;
end

% save all the variables as function of time in a text file for ploting later 
timeaxis = pop.x.';
save(strcat(filename, '-all_vs_t'),"timeaxis","*h","*m","Cl")

% plot all the variables as function of time and save them in a PDF file
plot_all(strcat(filename,'.pdf'),pop.x,Sh,Sp,Ih,Ip,Rh,Dh,Lm,Sm,Im,Cm,Cl,N);


% find the peak in the number of infected cases and save
Infmax = zeros(1,N);     % maximum number of infected people at the same time
peak_time = zeros(1,N);  % time at which the peak occurs
for jj=1:N
    [Infmax(jj), peak_time(jj)] = max( Ih(:,jj) );
    
end

% S_out and R_out only save the last value of each of each patch
% I_out saves the maximum value of infected individuals
% these variables are the return value of this function
S_out = Sh(end,:) ;
I_out = Infmax(:)';
R_out = Rh(end,:) ;

% calculate number of recovered peopple as a function 
% of the distance to the first infected patch (J0)
% dd(J0,i) is the distance from initial patch J0
% to any patch "i"
R_vs_d = zeros(1, lat_size);
cnt_patch = zeros(1, lat_size);
for i = 1:N
    if dd(J0,i) > 0 
        R_vs_d(dd(J0,i)) = R_vs_d(dd(J0,i)) + R_out(i);
        cnt_patch(dd(J0,i)) = cnt_patch(dd(J0,i)) + 1 ;
    end 
end


% save the time at which the peak of infected people occurs and 
% as a function of the distance to the 
delay_vs_d = [];
for i = 1:N
    delay_vs_d = [delay_vs_d; dd(J0,i), timeaxis(peak_time(i))];
end 

R_vs_d = R_vs_d ./ cnt_patch; % normalization

% save all those quantities on different text files 
dlmwrite(strcat(filename, '-STATIONARY-SUS.csv'), [tau, epsilon(1), mean(S_out)], 'delimiter',' ')
dlmwrite(strcat(filename, '-STATIONARY-INF.csv'), [tau, epsilon(1), mean(I_out)], 'delimiter',' ')
dlmwrite(strcat(filename, '-STATIONARY-REC.csv'), [tau, epsilon(1), mean(R_out)], 'delimiter',' ')
dlmwrite(strcat(filename, '-REC-vs-d.csv'), [[1:lat_size]', R_vs_d'], 'delimiter',' ')
dlmwrite(strcat(filename, '-peak_time-vs-d.csv'), delay_vs_d, 'delimiter',' ')


function dX = SIR_eq_delay_global(t, X, Z) % ** Using this one
    
    % This function uses all the infected in the system to 
    % calculate the amount of pesticides for each patch
    

    %m = 12; %% number of equations for each patch
    
    % References:
    % X(1+m*j) = fully susceptible humans (Sh)
    % X(2+m*j) = protected susceptible humans (Sp)
    % X(3+m*j) = infected humans (Ih)
    % X(4+m*j) = protected infected humans; (Ip)
    % X(5+m*j) = recovered humans (Rh)
    % X(6+m*j) = Severe Outcomes (Dh)
    % X(7+m*j) = mosquito larvae (Lm)
    % X(8+m*j) = susceptible mosquitoes (Sm)
    % X(9+m*j) = infected mosquitoes (Im)
    % X(10+m*j) = mosquito control (Cj)
    % X(11+m*j) = larvae control (Cl)
    % X(12+m*j) = total number of recovered people for 

    dX = zeros(m*N,1) ; % dX has to be a column vector

    
    for j = 0:N-1 
        
        N_total = X(1+m*j) + X(2+m*j) + X(3+m*j) + X(4+m*j) + X(5+m*j);
        
        Ih_lag = 0;
        Dh_lag = 0;
        for k = 0:N-1
            Ih_lag = Ih_lag + Z(3+m*k, k+1); 
            Dh_lag = Dh_lag + Z(6+m*k, 2*k+1);
        end
        
        flux = flux_of_people(X(1:m:end), j+1, mov, N);
        dX(1+j*m) = b*(X(1+m*j) + X(2+m*j) + X(3+m*j) + X(4+m*j) + X(5+m*j)) - beta(1)*X(9+m*j)*X(1+m*j) - (gamma_D*(X(3+m*j)+X(4+m*j))/N_total + gamma_B*(X(8+m*j)+X(9+m*j)))*X(1+m*j) - mu(1)*X(1+m*j) + lamda*X(2+m*j) + flux; %=dS^h/dt
        
        flux = flux_of_people(X(2:m:end), j+1, mov, N);
        dX(2+j*m) = (gamma_D*(X(3+m*j)+X(4+m*j))/N_total + gamma_B*(X(8+m*j)+X(9+m*j)))*X(1+m*j) - rho*beta(1)*X(9+m*j)*X(2+m*j) - mu(1)*X(2+m*j) - lamda*X(2+m*j) + flux; %dS^p/dt
        
        flux = flux_of_people(X(3:m:end), j+1, mov, N);
        %dX(3+j*m) = beta(1)*X(9+m*j)*X(1+m*j) - X(3+m*j)*(gamma_D*(X(3+m*j)+X(4+m*j))/N_total + gamma_B*(X(8+m*j)+X(9+m*j))) + X(3+m*j)*( -r - mu(1) - w ) + flux; %=dI^h/dt
        dX(3+j*m) = beta(1)*X(9+m*j)*X(1+m*j) + X(3+m*j)*( -r - mu(1) - w ) + flux; %=dI^h/dt
        
        flux = flux_of_people(X(4:m:end), j+1, mov, N);
        %dX(4+j*m) = (gamma_D*(X(3+m*j)+X(4+m*j))/N_total + gamma_B*(X(8+m*j)+X(9+m*j)))*X(3+m*j) - mu(1)*X(4+j*m) + rho*beta(1)*X(9+j*m)*X(2+j*m) - r*X(4+m*j) + flux;%=dI^p/dt        
        dX(4+j*m) = rho*beta(1)*X(9+j*m)*X(2+j*m) - mu(1)*X(4+j*m) - r*X(4+m*j) + flux;%=dI^p/dt        
        
        flux = flux_of_people(X(5:m:end), j+1, mov, N);
        dX(5+j*m) = r*X(3+m*j) + r*X(4+m*j) - mu(1)*X(5+m*j) + flux;%=dR^h/dt
        
        dX(6+j*m) = delta*X(3+m*j);%=dD^h/dt

        tmp = F(eta * (X(8+m*j) + X(9+m*j)), Kl(j+1)) - nu*X(7+m*j) - X(11+m*j)*X(7+m*j) ;  % = dL^m/dt    
        dX(7+j*m) = check_neg(X(7+m*j), tmp);
         
        tmp = nu*X(7+m*j) - beta(2)*X(3+m*j)*X(8+m*j) - rho*beta(2)*X(4+m*j)*X(8+m*j) - mu(2)*X(8+m*j) - X(8+m*j)*X(10+m*j); %=dS^m/dt
        dX(8+j*m) = check_neg(X(8+m*j), tmp);
        
        tmp = beta(2)*X(3+m*j)*X(8+m*j) + rho*beta(2)*X(4+m*j)*X(8+m*j) - mu(2)*X(9+m*j) - X(10+m*j)*X(9+m*j); %=dI^m/dt 
        dX(9+j*m) = check_neg(X(9+m*j), tmp);
        
        tmp = alpha(1)*Dh_lag + gamma(1)*Ih_lag - (epsilon(1)*X(10+m*j)); % mosquito control
        dX(10+j*m) = check_neg_or_one(X(10+m*j), tmp);
        
        tmp = alpha(2)*Dh_lag + gamma(2)*Ih_lag - (epsilon(2)*X(11+m*j)); % larvae control
        dX(11+j*m) = check_neg_or_one(X(11+m*j), tmp);
        
        dX(12+j*m) = r*X(3+m*j);

    end
end


function dX = SIR_eq_no_delay_global(t, X)

    m = 12; %% number of equations for each patch
    
    % References:
    % X(1+m*j) = susceptible humans
    % X(2+m*j) = infected humans
    % X(3+m*j) = recovered humans
    % X(4+m*j) = Severe Outcomes
    % X(5+m*j) = mosquito larvae
    % X(6+m*j) = susceptible mosquitoes
    % X(7+m*j) = infected mosquitoes
    % X(8+m*j) = mosquito control
    % X(9+m*j) = larvae control
    % X(10+m*j) = total number of recovered people for 
    
    %X = pop(1:m*N);
    dX = zeros(m*N,1);
    
    for j = 0:N-1 
        
        Ih_tot = sum(X(3:m:end));
        Dh_tot = sum(X(6:m:end));
        
        flux = flux_of_people(X(1:m:end), j+1, mov, N);
        dX(1+j*m) = b*(X(1+m*j) + X(2+m*j) + X(3+m*j) + X(4+m*j) + X(5+m*j)) - beta(1)*X(9+m*j)*X(1+m*j) - (gamma_D*(X(3+m*j)+X(4+m*j)) + gamma_B*(X(8+m*j)+X(9+m*j)))*X(1+m*j) - mu(1)*X(1+m*j) + lamda*X(2+m*j) + flux; %=dS^h/dt
        
        flux = flux_of_people(X(2:m:end), j+1, mov, N);
        dX(2+j*m) = (gamma_D*(X(3+m*j)+X(4+m*j)) + gamma_B*(X(8+m*j)+X(9+m*j)))*X(1+m*j) - rho*beta(1)*X(9+m*j)*X(2+m*j) - mu(1)*X(2+m*j) - lamda*X(2+m*j) + flux; %dS^p/dt
        
        flux = flux_of_people(X(3:m:end), j+1, mov, N);
        dX(3+j*m) = beta(1)*X(9+m*j)*X(1+m*j) + lamda*X(4+m*j) + X(3+m*j)*( -r - mu(1) - w ) + flux; %=dI^h/dt
        
        flux = flux_of_people(X(4:m:end), j+1, mov, N);
        dX(4+j*m) = (gamma_D*(X(3+m*j)+X(4+m*j)) + gamma_B*(X(8+m*j)+X(9+m*j)))*X(3+m*j) - lamda*X(4+j*m) - mu(1)*X(4+j*m) + rho*beta(1)*X(9+j*m)*X(2+j*m) + flux;%=dI^p/dt        
        
        flux = flux_of_people(X(5:m:end), j+1, mov, N);
        dX(5+j*m) = r*X(3+m*j) - mu(1)*X(5+m*j) + flux; %=dR^h/dt
        
        dX(6+j*m) = delta*X(3+m*j);%=dD^h/dt

        tmp = F(eta * (X(8+m*j) + X(9+m*j)), Kl(j+1)) - nu*X(7+m*j) - X(11+m*j)*X(7+m*j) ;  % = dL^m/dt    
        dX(7+j*m) = check_neg(X(7+m*j), tmp);
         
        tmp = nu*X(7+m*j) - beta(2)*X(3+m*j)*X(8+m*j) - rho*beta(2)*X(4+m*j)*X(8+m*j) - mu(2)*X(8+m*j) - X(8+m*j)*X(10+m*j); %=dS^m/dt
        dX(8+j*m) = check_neg(X(8+m*j), tmp);
        
        tmp = beta(2)*X(3+m*j)*X(8+m*j) + rho*beta(2)*X(4+m*j)*X(8+m*j) - mu(2)*X(9+m*j) - X(10+m*j)*X(9+m*j); %=dI^m/dt 
        dX(9+j*m) = check_neg(X(9+m*j), tmp);
        
        tmp = alpha(1)*Dh_tot + gamma(1)*Ih_tot  - (epsilon(1)*X(10+m*j)); % mosquito control
        dX(10+j*m) = check_neg_or_one(X(10+m*j), tmp);
        
        tmp = alpha(2)*Dh_tot + gamma(2)*Ih_tot - (epsilon(2)*X(11+m*j)); % larvae control
        dX(11+j*m) = check_neg_or_one(X(11+m*j), tmp);
        
        dX(12+j*m) = r*X(3+m*j);

    end
    
end

end



function adj = square_adj_mat(N)
% calculate the adjacency matrix for a square lattice with connections
% to the first nighbohrs only, and fixed boundary conditions
% Set the matrix size
r = N;  % number of rows
c = N;  % number of columns                            
diagVec1 = repmat([ones(c-1, 1); 0], r, 1);  % Make the first diagonal vector
                                             %   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);                % Remove the last value
diagVec2 = ones(c*(r-1), 1);                 % Make the second diagonal vector
                                             %   (for vertical connections)
adj = diag(diagVec1, 1)+diag(diagVec2, c);   % Add the diagonals to a zero matrix
adj = sparse(adj+adj.');   

end

function new = check_neg_or_one(old_value, change)
% check if a variable is about to be larger than one

    if  old_value + change > 1
        new = 1 - old_value;
        return 
    else
        new = change;
        return 
    end
% check if a variable is about to be lower than zero
    if  old_value + change < 0
        new = -old_value;
        return 
    else
        new = change;
        return 
    end
end 

function new = check_neg(old_value, change)
% check if a variable is about to be lower than zero
    if  old_value + change < 0
        new = -old_value;
    else
        new = change;
    end
end 

function cc = F(x, Kl)
% function for carrying capacity in the eq. for larvaes:
    cc = x * ( 1 - ( x / Kl ));
    %if cc < 0
        %cc = 0;  
    %end
end

function flux = flux_of_people(X, current_patch, mov, N)
% calculate the term for the in and out flux of people in a patch
% mov(i,j) is the fraction of people that moves from i to j

    j = current_patch;

    in_flux_tmp = mov * X;
    in_flux = in_flux_tmp(j);
    
    out_flux = sum(mov(j,:)) * X(j) ;
    
    flux = in_flux - out_flux ;
end

function h = plot_all(txt,t,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,N)
fig=figure('Name','SIR', 'units','inch','position',[0,0,14,8]) ;

fig.PaperType = 'uslegal';
fig.PaperOrientation = 'landscape';

subplot(3,3,1);
h1=plot(t,X1(:,1),'r');
hold on
h=plot(t,X1(:,1:N),'r');
h2=plot(t,X2(:,1),'b');
h=plot(t,X2(:,1:N),'b');
xlabel 'Time';
ylabel 'Susceptible Humans';
legend([h1, h2],'Unprotected','Protected')

subplot(3,3,2);
h1=plot(t,X3(:,1),'r');
hold on
h=plot(t,X3(:,1:N),'r');
h=plot(t,X4(:,1:N),'b');
h2=plot(t,X4(:,1),'b');
xlabel 'Time';
ylabel 'Infected Humans';
%legend([h1, h2],'Non-protected Infected','Protected Infected','Location','northeast')

title(txt)

subplot(3,3,3);
h=plot(t,X5(:,1:N));
xlabel 'Time';
ylabel 'Recovered Humans';

subplot(3,3,4) ;
h=plot(t,X6(:,1:N));
xlabel 'Time';
ylabel 'Severe Outcomes';

subplot(3,3,5) ;
h=plot(t,X10(:,1:N));
xlabel 'Time';
ylabel 'Mosquitoes Control';

subplot(3,3,6) ;
h=plot(t,X11(:,1:N));
xlabel 'Time';
ylabel 'Larvae Control';

subplot(3,3,7) ;
h=plot(t,X7(:,1:N));
xlabel 'Time';
ylabel 'Larvae';

subplot(3,3,8) ;
h=plot(t,X8(:,1:N));
xlabel 'Time';
ylabel 'Susceptible Mosquitoes';

subplot(3,3,9) ;
h=plot(t,X9(:,1:N));
xlabel 'Time';
ylabel 'Infected Mosquitoes';

set(gcf,'color','white')

%print(txt,'-depsc')
saveas(gcf, txt);  %txt is the path and filename 

end

function D=distance(A)
% calculate de distance matrix from a given adjancency matrix A

N = max(size(A));
D = NaN(N);
B = A;
k = 1;
%while any(isnan(D(:)))
while k <= N
    % Check for new walks, and assign distance
    D(B>0 & isnan(D)) = k;
    % Iteration
    k = k + 1;
    B = B * A;
end
D(1:N+1:end) = 0; % diagonal elements are always 0
% Now D contains the distance matrix
end 
