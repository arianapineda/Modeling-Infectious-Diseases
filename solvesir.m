% Ariana Pineda, CAAM 210, SPRING 2022, Modeling Infectious Diseases
% solvesir.m
% this script models infectious diseases in a population during a given
% time frame using information on infectious and recovered people. 
% Last modified: February 16, 2022

function solvesir
% this driver initializes parameters for each SIR model and generates plots
% for each model.

% part 1: initialize parameters
alpha = 0.7;
beta = 0.1;
M = 7.9e6;
R0 = 0;
I0 = 10;
initialval = [M-I0-R0, R0, I0];
Tfinal = 150;
t = 0:150;
[Sval,Rval,Ival] = simpleSIR(M, alpha, beta, initialval, Tfinal);

% part 1: plot
figure(1)
hold on
plot(t, Sval, 'LineWidth', 2)
plot(t, Rval, 'LineWidth', 2)
plot(t, Ival, 'LineWidth', 2)
legend('Susceptibles','Recovered','Infectious')
xlabel('nb of days')
ylabel('population')
hold off

% part 2: initialize parameters
alpha = 0.5;
beta = 1/3;
gamma = 0.01;
mu = 1/(76*365);
S0 = 7.9e6;
R0 = 0;
I0 = 10;
Tfinal = 4*365;
t = 0:Tfinal;
initialval = [S0, R0, I0];

[Sval,Rval,Ival,M] = variableSIR(alpha,beta,gamma,mu,initialval,Tfinal);

% part 2: plot 1
figure(2)
hold on
plot(t, Sval, 'LineWidth', 2)
plot(t, Rval, 'LineWidth', 2)
plot(t, Ival, 'LineWidth', 2)
plot(t,M,'LineWidth',2)
xlabel('nb of days')
ylabel('population')
legend('Susceptibles', 'Recovered','Infectious','Total')
hold off

% part 2: plot 2
figure(3)
hold on
plot(Sval,Ival,'LineWidth',2)
xlabel('Susceptibles')
ylabel('Infected')
hold off

% part 3: initialize parameters

omega = 1/365;
Tfinal = 4*365;
t = 0:Tfinal;
[Sval,Rval,Ival,M] = variableimmSIR(alpha,beta,gamma,mu,initialval,Tfinal, omega);


% part 3: plot 1
figure(4)
hold on
plot(t, Sval, 'LineWidth', 2)
plot(t, Rval, 'LineWidth', 2)
plot(t, Ival, 'LineWidth', 2)
plot(t,M,'LineWidth',2)
xlabel('nb of days')
ylabel('population')
legend('Susceptibles', 'Recovered','Infectious','Total')
hold off

% part 3: plot 2
figure(5)
hold on
plot(Sval,Ival,'LineWidth', 2)
xlabel('Susceptibles')
ylabel('Infected')
hold off
end

function [Sval,Rval,Ival] = simpleSIR(M,alpha,beta,initialval,Tfinal)
% creates a SIR model using given parameters over a period of time

% inputs: M-total population, alpha - contacts per day sufficient to spread
% disease, beta - fraction of infected group that will recover during any 
% given day, initialval - vector of initial values of S, I, and R, 
% Tfinal - time in days

%outputs: Sval-number of susceptible individuals; these are healthy
%individuals, Rval-number of recovered individuals, Ival-number of infected
%individuals


% initialization of the variables
Sval(1) = initialval(1);
Rval(1) = initialval(2);
Ival(1) = M-Sval(1)-Rval(1);

%iterates over specified time period
for i=1:Tfinal
    Sval(i+1)=Sval(i)-((alpha/M)*Sval(i)*Ival(i));
    Rval(i+1)=Rval(i)+(beta*Ival(i));
    Ival(i+1) = M-Sval(i+1)-Rval(i+1);
end
end


%Interpret the results from Fig(1) and your new fig-
%ure, i.e. why S(t),I(t),R(t) changed in a specific way when we changed
%the parameters to receive the new figure.

%In the new figure, the susceptibles drastically decreased as alpha increased and beta increased, the recovered drastically increased,
%and the number of infectious spiked higher and earlier. That means
%if there are more contacts per day suffient to spread the disease, people get infected 
%more quickly. With a higher beta, more people will recover. The reason for this is intuitive.
%If more people are infected, more people can spread the disease quicker. This will lead to a quick drop in
%susceptible people since they are probably infectious. As the number of infected people 
%increase and approaches the whole population, the number of recoveries will rise drastically because there are
%more people who have the ability to recover. You can't recover from a disease if you don't have one.
%As more people recover, the number of infected drops.

function [Sval, Rval, Ival, M] = variableSIR(alpha,beta,gamma,mu,initialval,Tfinal)
% creates a variable SIR model using given parameters

% inputs: alpha - contacts per day sufficient to spread disease, beta - 
% fraction of infected group that will recover during any given day, gamma 
% - rate of loss of individuals, mu - death rate for all population from 
% other causes, initialval - vector of initial values of S, I, and R, 
% Tfinal - time in days

% outputs: Sval-number of susceptible individuals; these are healthy
%individuals, Rval-number of recovered individuals, Ival-number of infected
%individuals, M-total population

%initialize variables
Sval(1) = initialval(1);
Rval(1) = initialval(2);
Ival(1) = initialval(3);
M(1) = Sval(1)+Rval(1)+Ival(1);
% loop over the time steps
for i=1:Tfinal
    Sval(i+1)=Sval(i)-((alpha/M(i))*Sval(i)*Ival(i))+(mu*M(i))-(mu*Sval(i));
    Rval(i+1)=Rval(i)+(beta*Ival(i))-(mu*Rval(i));
    Ival(i+1) = Ival(i)+((alpha/M(i))*Sval(i)*Ival(i))-(beta*Ival(i))-(mu+gamma)*Ival(i);
    M(i+1) = Sval(i+1) + Rval(i+1) + Ival(i+1);
end
end

function [Sval, Rval, Ival, M] = variableimmSIR(alpha,beta,gamma,mu,initialval,Tfinal,omega)
% creates a modified variable SIR model using given parameters that takes into
% account loss of immunity

% inputs: alpha - contacts per day sufficient to spread disease, beta - 
% fraction of infected group that will recover during any given day, gamma 
% - rate of loss of individuals, mu - death rate for all population from 
% other causes, initialval - vector of initial values of S, I, and R, 
% Tfinal - time in days, omega - a fraction of the recovered population 
% moves to the group of susceptibles

% outputs: Sval-number of susceptible individuals; these are healthy
%individuals, Rval-number of recovered individuals, Ival-number of infected
%individuals, M-total population

%initialize variables
Sval(1) = initialval(1);
Rval(1) = initialval(2);
Ival(1) = initialval(3);
M(1) = Sval(1) + Rval(1) + Ival(1);

%iterates over specified time period
for i = 1:Tfinal
    Sval(i+1)=Sval(i)-((alpha/M(i))*Sval(i)*Ival(i))+(mu*M(i))-(mu*Sval(i))+omega*Rval(i);
    Rval(i+1)=Rval(i)+(beta*Ival(i))-(mu*Rval(i))-omega*Rval(i);
    Ival(i+1) = Ival(i)+((alpha/M(i))*Sval(i)*Ival(i))-(beta*Ival(i))-(mu+gamma)*Ival(i);
    M(i+1) = Sval(i+1) + Rval(i+1) + Ival(i+1);
end 
end