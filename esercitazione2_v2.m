clear all;
clc;

%% Esercizio n*1 

% Dati modello 
V2 = 5;
k01 = 1.2;
k02 = 1.2;
Vmax = 110;
km = 50;
Dv = 500;

% Risoluzione sistema eq. differenziali
tspan = [0 10];
q0 = [Dv 0];
[t,q] = ode45(@(t,q) odefcn(t,q,k01,k02,Vmax,km), tspan, q0);

% Plot delle quantita
figure(1);
subplot(1,2,1), plot(t,q(:,1)), title("Quantita compartimento 1"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on;
subplot(1,2,2), plot(t,q(:,2)), title("Quantita compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on;

% Plot della concentrazione 
figure(2);
c2 = q(:,2)/V2;
plot(t,c2), title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on;

% Plot flusso
figure(3);
vel = [];
for i = 1:length(q)
    vel(i) = (Vmax/(km+q(i,1)))*q(i,1);
end
plot(t,vel), title("Grafico velocita di assorbimento"),
xlabel("Tempo (ore)"), ylabel("Velocita (mg/ora)"), grid on;

%% Punto A - Confronto con sistema lineare

% Confronto con modello lineare
% Dati modello 
kl01 = 1.2;
kl02 = 1.2;
kl21 = 2.2;
V2 = 5;
d = 500;

% Matrici del sistema compartimentale lineare
A = [-(kl01+kl21),0;kl21,-kl02]; 
B = [d;0];
C = [0,1/V2];
D = 0;

% Creazione del modello compartimentale lineare
sys = ss(A,B,C,D);

figure(4);
t1 = [0:0.01:10];
[ct,T,qt] = impulse(sys,t1);
subplot(1,2,1), plot(T,qt), hold on, plot(t,q(:,1),t,q(:,2)), hold off, title("Quantita nel compartimento"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), 
legend("Compartimento 1 - Lineare","Compartimento 2 - Lineare","Compartimento 1 - MM","Compartimento 2 - MM"), grid on;
subplot(1,2,2), plot(T,ct), hold on, plot(t,c2), hold off, title("Concentrazione nel compartimento 2"),
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), legend("Concentrazione - Lineare", "Concentrazione - MM"), grid on;

%% Punto B - Variazione Vmax e km 
% Varizione Vmax
Vmax = Vmax*2;

% Risoluzione sistema eq. differenziali
tspan = [0 10];
q0 = [Dv 0];
[t_v,q_v] = ode45(@(t,q) odefcn(t,q,k01,k02,Vmax,km), tspan, q0);

% Plot delle quantita
figure(5);
subplot(1,2,1), plot(t,q(:,1)), hold on, plot(t_v,q_v(:,1)), hold off, 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 1 - Vmax = 110","Compartimento 1 - Vmax = 220"),
title("Quantita compartimento 1");
subplot(1,2,2), plot(t,q(:,2)), hold on, plot(t_v,q_v(:,2)), hold off,
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 2 - Vmax = 110","Compartimento 2 - Vmax = 220"),
title("Quantita compartimento 2");

% Plot della concentrazione 
figure(6);
c2 = q(:,2)/V2;
c2_v = q_v(:,2)/V2;
plot(t,c2), hold on, plot(t_v,c2_v), hold off, title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on,
legend("Concentrazione - Vmax = 110","Concentrazione - Vmax = 220");

% Variazione km
km = km/2;

% Risoluzione sistema eq. differenziali
tspan = [0 10];
q0 = [Dv 0];
[t_k,q_k] = ode45(@(t,q) odefcn(t,q,k01,k02,Vmax,km), tspan, q0);

% Plot delle quantita
figure(7);
subplot(1,2,1), plot(t,q(:,1)), hold on, plot(t_k,q_k(:,1)), hold off, 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 1 - km = 50","Compartimento 1 - km = 25"),
title("Quantita compartimento 1");
subplot(1,2,2), plot(t,q(:,2)), hold on, plot(t_k,q_k(:,2)), hold off,
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 2 - km = 50","Compartimento 2 - km = 25"),
title("Quantita compartimento 2");

% Plot della concentrazione 
figure(8);
c2 = q(:,2)/V2;
c2_k = q_k(:,2)/V2;
plot(t,c2), hold on, plot(t_k,c2_k), hold off, title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on,
legend("Concentrazione - km = 50","Concentrazione - km = 25");

%% Punto C - Cambiamento dose somministrata
% Dati modello 
V2 = 5;
k01 = 1.2;
k02 = 1.2;
Vmax = 110;
km = 50;
Dv = 500;

% Risoluzione sistema eq. differenziali
tspan = [0,10];
q0 = [(Dv/10),0];
[t,q] = ode45(@(t,q) odefcn(t,q,k01,k02,Vmax,km), tspan, q0);

c2 = q(:,2)/V2;

figure(9);
t1 = [0:0.01:10];
[ct,T,qt] = impulse(sys,t1);
subplot(1,2,1), plot(T,qt), hold on, plot(t,q(:,1),t,q(:,2)), hold off, title("Quantita nel compartimento - Dv = 50"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), 
legend("Compartimento 1 - Lineare","Compartimento 2 - Lineare","Compartimento 1 - MM","Compartimento 2 - MM"), grid on;
subplot(1,2,2), plot(T,ct), hold on, plot(t,c2), hold off, title("Concentrazione nel compartimento 2 - Dv = 50"),
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), legend("Concentrazione - Lineare", "Concentrazione - MM"), grid on;

%% Punto D - Cinetica Hill
% Dati modello 
V2 = 5;
k01 = 1.2;
k02 = 1.2;
Vmax = 110;
km = 50;
Dv = 500;

% Vettore esponenti della cinetica di Hill
e = [1:0.5:3];

tspan = [0,10];
q0 = [Dv,0];

figure(10);

for i = 1:length(e)
    esp = ['Hill q = ' num2str(e(i))];
    [t,q] = ode45(@(t,q) odefcnHill(t,q,e(i),k01,k02,Vmax,km), tspan, q0);
    subplot(1,2,1), hold on, plot(t,q(:,1),'DisplayName',esp), title("Quantita nel compartimento 1"), 
    xlabel("Tempo (ore)"), ylabel("Quantita (mg/L)"), grid on, legend("show");
    subplot(1,2,2), hold on, plot(t,q(:,2),'DisplayName',esp), title("Quantita nel compartimento 2"), 
    xlabel("Tempo (ore)"), ylabel("Quantita (mg/L)"), grid on, legend("show");
end

figure(11);

for i = 1:length(e)
    esp = ['Hill q = ' num2str(e(i))];
    [t,q] = ode45(@(t,q) odefcnHill(t,q,e(i),k01,k02,Vmax,km), tspan, q0);
    c = zeros(length(q));
    c = q(:,2)/V2;
    hold on, plot(t,c,'DisplayName',esp), title("Concentrazione compartimento 2"), 
    xlabel("Tempo (ore)"), ylabel("Quantita (mg/L)"), grid on, legend("show");
end



