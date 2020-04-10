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

% Q.ta nel primo compartimento
tspan1 = [0,10];
q1_0 = Dv;
[t1,q1] = ode45(@(t1,q1) -(k01+Vmax/(km+q1))*q1,tspan1,q1_0);

% Q.ta nel secondo compartimento
tspan2 = [0,10];
q2_0 = 0;
[t2,q2] = ode45(@(t2,q2) (Vmax/(km+q1))*q1-k02*q2,tspan2,q2_0);

% Plot delle quantita
figure(1);
subplot(2,1,1), plot(t1,q1), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on;
subplot(2,1,2), plot(t2,q2), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on;

% Plot della concentrazione 
figure(2);
c2 = q2/V2;
plot(t2,c2), title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on;

% Plot della velocita
figure(3);
vel = [];
for i = 1:length(q1)
    vel(i) = (Vmax/(km+q1(i)))*q1(i);
end
plot(t1,vel), title("Grafico velocita di assorbimento"),
xlabel("Tempo (ore)"), ylabel("Velocita (mg/ora)"), grid on;

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

%% Punto A - Grafici della quantita e della concentrazione
figure(4);
t = [0:0.01:10];
[ct,T,qt] = impulse(sys,t);
subplot(1,2,1), plot(T,qt), hold on, plot(t1,q1,t2,q2), hold off, title("Quantita nel compartimento"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), 
legend("Compartimento 1 - Lineare","Compartimento 2 - Lineare","Compartimento 1 - MM","Compartimento 2 - MM"), grid on;
subplot(1,2,2), plot(T,ct), hold on, plot(t2,c2), hold off, title("Concentrazione nel compartimento 2"),
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), legend("Concentrazione - Lineare", "Concentrazione - MM"), grid on;


%% Punto B - Variazione Vmax e km 
% Varizione Vmax
Vmax = Vmax*2;

% Q.ta nel primo compartimento
tspan1_2 = [0,10];
q1_0_2 = Dv;
[t1_2,q1_2] = ode45(@(t1_2,q1_2) -(k01+Vmax/(km+q1_2))*q1_2,tspan1_2,q1_0_2);

% Q.ta nel secondo compartimento
tspan2_2 = [0,10];
q2_0_2 = 0;
[t2_2,q2_2] = ode45(@(t2_2,q2_2) (Vmax/(km+q1_2))*q1_2-k02*q2_2,tspan2_2,q2_0);

% Plot delle quantita
figure(5);
subplot(1,2,1), plot(t1,q1), hold on, plot(t1_2,q1_2), hold off, 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 1 - Vmax = 110","Compartimento 1 - Vmax = 220");
subplot(1,2,2), plot(t2,q2), hold on, plot(t2_2,q2_2), hold off,
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 2 - Vmax = 110","Compartimento 2 - Vmax = 220");

% Plot della concentrazione 
figure(6);
c2 = q2/V2;
c2_2 = q2_2/V2;
plot(t2,c2), hold on, plot(t2_2,c2_2), hold off, title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on,
legend("Concentrazione - Vmax = 110","Concentrazione - Vmax = 220");

% Variazione km
km = km/2;

% Q.ta nel primo compartimento
tspan1_2 = [0,10];
q1_0_2 = Dv;
[t1_2,q1_2] = ode45(@(t1_2,q1_2) -(k01+Vmax/(km+q1_2))*q1_2,tspan1_2,q1_0_2);

% Q.ta nel secondo compartimento
tspan2_2 = [0,10];
q2_0_2 = 0;
[t2_2,q2_2] = ode45(@(t2_2,q2_2) (Vmax/(km+q1_2))*q1_2-k02*q2_2,tspan2_2,q2_0);

% Plot delle quantita
figure(7);
subplot(1,2,1), plot(t1,q1), hold on, plot(t1_2,q1_2), hold off, 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 1 - km = 50","Compartimento 1 - km = 25");
subplot(1,2,2), plot(t2,q2), hold on, plot(t2_2,q2_2), hold off,
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on,
legend("Compartimento 2 - km = 50","Compartimento 2 - km = 25");

% Plot della concentrazione 
figure(8);
c2 = q2/V2;
c2_2 = q2_2/V2;
plot(t2,c2), hold on, plot(t2_2,c2_2), hold off, title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on,
legend("Concentrazione - km = 50","Concentrazione - km = 25");


%% Punto C - Cambiamento dose somministrato
% Q.ta nel primo compartimento
tspan1 = [0,10];
q1_0 = Dv/10;
[t1,q1] = ode45(@(t1,q1) -(k01+Vmax/(km+q1))*q1,tspan1,q1_0);

% Q.ta nel secondo compartimento
tspan2 = [0,10];
q2_0 = 0;
[t2,q2] = ode45(@(t2,q2) (Vmax/(km+q1))*q1-k02*q2,tspan2,q2_0);

% Plot della concentrazione 
figure(2);
c2 = q2/V2;
plot(t2,c2), title("Concentrazione compartimento 2"), 
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on;

figure(9);
t = [0:0.01:10];
[ct,T,qt] = impulse(sys,t);
subplot(1,2,1), plot(T,qt), hold on, plot(t1,q1,t2,q2), hold off, title("Quantita nel compartimento - Dv = 50"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), 
legend("Compartimento 1 - Lineare","Compartimento 2 - Lineare","Compartimento 1 - MM","Compartimento 2 - MM"), grid on;
subplot(1,2,2), plot(T,ct), hold on, plot(t2,c2), hold off, title("Concentrazione nel compartimento 2 - Dv = 50"),
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), legend("Concentrazione - Lineare", "Concentrazione - MM"), grid on;








