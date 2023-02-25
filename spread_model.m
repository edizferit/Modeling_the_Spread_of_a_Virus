% ME303 TERM PROJECT QUESTION 1
clear all, close all, clc

% 1: MODELING THE SPREAD OF A VIRUS

%============= GIVENS =============%
% STEP SIZE
h = 1/500;
% INTERVAL
tmin = 0;
tmax = 100; % days
t = tmin:h:tmax;
% NUMBER OF ITERATIONS
m = (tmax-tmin)/h;

%=========== PARAMETERS ===========%
c = 4; % encounters per day
B = 0.2; % transmission probability per encounter
a = 0.125; % rate at which infected become infectious (per day)
g = 0.1; % rate at which an infected person becomes symptomatic (per day)
w = 0.2; % recovery rate
AB = 1000; % HOSPITAL BEDS

%========== SIMPLIFIED MATHEMATICAL MODEL ==========%

% 5 SUB GROUPS OF A POPULATION
S = zeros(1,m); % Susceptibles
E = zeros(1,m); % Exposed
I = zeros(1,m); % Infected
M = zeros(1,m); % Medically Symptomatic
R = zeros(1,m); % Recovered

% INITIAL CONDITIONS
S(1) = 10000;
E(1) = 10;
I(1) = 0;
M(1) = 0;
R(1) = 0;

% DIFFERENTIAL EQUATIONS
dS = @(S,E,I,M,R) -c*B*(I./(S+E+I+M+R)).*S; % S prime
dE = @(S,E,I,M,R) c*B*(I./(S+E+I+M+R)).*S - a*E; % E prime
dI = @(S,E,I,M,R) a*E - g*I; % I prime
dM = @(S,E,I,M,R) g*I - w*M; % M prime
dR = @(S,E,I,M,R) w*M; % R prime

%=========== MORE COMPACT FORM OF THE MODEL ===========%

% function handle
dy{1} = @(V) dS(V(1),V(2),V(3),V(4),V(5));
dy{2} = @(V) dE(V(1),V(2),V(3),V(4),V(5));
dy{3} = @(V) dI(V(1),V(2),V(3),V(4),V(5));
dy{4} = @(V) dM(V(1),V(2),V(3),V(4),V(5));
dy{5} = @(V) dR(V(1),V(2),V(3),V(4),V(5));

% vector matrix
V = [S; E; I; M; R];


% ======== EULER SOLUTION ========%
for i=1:m
    for j=1:5
        V(j,i+1) = V(j,i) + h*dy{j}(V(:,i));
    end
end

S_e = V(1,:); % Store the euler values
E_e = V(2,:);
I_e = V(3,:);
M_e = V(4,:);
R_e = V(5,:);
V_e = [S_e; E_e; I_e; M_e; R_e];

figure('Position', [50 100 1400 550]) % create the common figure
subplot(1,3,1)
plot(t,S_e, 'r', 'linewidth', 2)
hold on, grid on
plot(t,E_e, 'g', 'linewidth', 2)
plot(t,I_e, 'b', 'linewidth', 2)
plot(t,M_e, 'y', 'linewidth', 2)
plot(t,R_e, 'k', 'linewidth', 2)
plot([tmin, tmax],[AB, AB], '--m', 'linewidth', 2)
title('Euler Solution')
xlabel('Time (days)')
ylabel('Population')
h2=legend('Susceptible (S)', 'Exposed (E)', 'Infected (I)',...
    'Medically Symptomatic (M)', 'Recovered (R)',...
    'Available Beds (AB)');
set(h2, 'Location', 'southoutside')


%======== 4TH ORDER RUNGE KUTTA SOLUTION ========%
f{1} = @(g,V,h,fv) g(V);
f{2} = @(g,V,h,fv) g(V + h/2*fv);
f{3} = @(g,V,h,fv) g(V + h/2*fv);
f{4} = @(g,V,h,fv) g(V + h*fv);
fv = zeros(5,4); % Each row stores f values for each equation
for i=1:m
    for j=1:4 % four f values
        for k=1:5 % five equations
            if j>1
                fv(k,j)= f{j}(dy{k}, V(:,i), h, fv(:,j-1)); % for f2, f3, f4
            else
                fv(k,j)= f{j}(dy{k}, V(:,i), h, fv(:,j)); % for f1
            end
        end
    end
    V(:,i+1) = V(:,i) + h/6*(fv(:,1) + 2*fv(:,2) + 2*fv(:,3) + fv(:,4));
end

S_rk = V(1,:); % Store the rk4 values
E_rk = V(2,:);
I_rk = V(3,:);
M_rk = V(4,:);
R_rk = V(5,:);
V_rk = [S_rk; E_rk; I_rk; M_rk; R_rk];

subplot(1,3,2)
plot(t,S_rk, 'r', 'linewidth', 2)
hold on, grid on
plot(t,E_rk, 'g', 'linewidth', 2)
plot(t,I_rk, 'b', 'linewidth', 2)
plot(t,M_rk, 'y', 'linewidth', 2)
plot(t,R_rk, 'k', 'linewidth', 2)
plot([tmin, tmax],[AB, AB], '--m', 'linewidth', 2)
title('Runge-Kutta 4 Solution')
xlabel('Time (days)')
ylabel('Population')
h1=legend('Susceptible (S)', 'Exposed (E)', 'Infected (I)',...
    'Medically Symptomatic (M)', 'Recovered (R)',...
    'Available Beds (AB)');
set(h1, 'Location', 'southoutside')


%======== ODE45 ========%
tspan = [tmin tmax];
y0 = [10000, 10, 0, 0, 0];
[t,y] = ode45(@(t,y) odefcn(y,c,B,a,g,w), tspan, y0);

S_ode = y(:,1)'; % Store the ode45 values (in row vectors)
E_ode = y(:,2)';
I_ode = y(:,3)';
M_ode = y(:,4)';
R_ode = y(:,5)';
V_ode = [S_ode; E_ode; I_ode; M_ode; R_ode];

subplot(1,3,3)
plot(t,y(:,1),'r',t,y(:,2),'g',t,y(:,3),'b',t,y(:,4),'y',...
    t,y(:,5),'k',[tmin tmax],[AB AB],'--m','linewidth',2)
grid on
title('ODE45 Solution')
xlabel('Time (days)')
ylabel('Population')
h3=legend('Susceptible (S)', 'Exposed (E)', 'Infected (I)',...
    'Medically Symptomatic (M)', 'Recovered (R)',...
    'Available Beds (AB)');
set(h3, 'Location', 'southoutside')


% The function that represents the system of equations for ode45
function dydt = odefcn(y,c,B,a,g,w)
dydt = zeros(5,1);
dydt(1) = -c*B*(y(3)./(y(1)+y(2)+y(3)+y(4)+y(5))).*y(1);
dydt(2) = c*B*(y(3)./(y(1)+y(2)+y(3)+y(4)+y(5))).*y(1) - a*y(2);
dydt(3) = a*y(2) - g*y(3);
dydt(4) = g*y(3) - w*y(4);
dydt(5) = w*y(4);
end
