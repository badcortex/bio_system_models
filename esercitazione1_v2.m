clear all;
clc;

%% Esercizio n*1 - Modello monocompartimentale

% Dati del modello 
ke = 1.2;
d = 500;
V1 = 5;

% Matrici del sistema compartimentale lineare
A = -ke;
B = d;
C = 1/V1;
D = 0;

% Creazione del modello compartimentale lineare
sys = ss(A,B,C,D);

% Grafici della quantita e della concentrazione
figure(1);
t = [0:0.01:10];
[ct,T,qt] = impulse(sys,t);
t_auc = [0, 0.76, 1.16, 1.72, 2.53, 4.01, 6];
c_auc = [100, 40.176, 24.8578, 12.6945, 4.8027, 0.81316, 0.074];
subplot(1,2,1),plot(T,qt), title("Quantita nel compartimento"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on;
subplot(1,2,2), plot(T,ct), hold on, plot(t_auc,c_auc,'o'), hold off, title("Concentrazione nel compartimento"),
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on;

% Grafico della concentrazione in scale logaritmica
figure(2);
semilogy(T,ct), title("Concentrazione compartimento (scala log)"),
xlabel("Tempo"), ylabel("Concentrazione"), grid on;

% Il sistema e' monocompartimentale e quindi ad una sola costante di tempo
% (e' possibile verificarlo dal grafico su scala logaritmica). Il tempo di 
% eliminazione e' pari a 5tau tau = 1/ke. 
tau = 1/ke;
te = 5*tau;

% Calcolo emivita 
tm = 0.693/ke;

% Poiche' il modello e' monocopartimentale concentrazione massima e
% concentrazione iniziale coincidono. 
c0 = d/V1;
cmax = d/V1;

% AUC - Calcolata con formula esatta e tramite approssimazione
auc1_1 = c0/ke;
auc1_2 = trapz(T,ct);

% MRT - Tempo medio di residenza 
% Uso l'approccio non compartimentale calcolando AUMC e AUC con
% campionamento 
auc1_nc = trapz(t_auc,c_auc);
aumc1_nc = trapz(t_auc,t_auc.*c_auc);

mrt = aumc1_nc/auc1_nc;
 
% CLtot - Clearance totale 
CLtot = V1*ke;


%% Esercizio n*2 - Modello con assorbimento 

% Dati modello 
k01 = 1.2;
k02 = 1.2;
k21 = 2.2;
V2 = 5;
d = 500;

% Matrici del sistema compartimentale lineare
A = [-(k01+k21),0;k21,-k02]; 
B = [d;0];
C = [0,1/V2];
D = 0;

% Creazione del modello compartimentale lineare
sys = ss(A,B,C,D);

% Grafici della quantita e della concentrazione
figure(3);
t = [0:0.01:10];
[ct,T,qt] = impulse(sys,t);
subplot(1,2,1), plot(T,qt), title("Quantita nel compartimento"), 
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), legend("Compartimento 1","Compartimento 2"), grid on;
subplot(1,2,2), plot(T,ct), title("Concentrazione nel compartimento 2"),
xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on;

% Grafico della concentrazione in scale logaritmica
figure(4);
semilogy(T,ct), title("Concentrazione 2* compartimento (scala log)"),
xlabel("Tempo"), ylabel("Concentrazione"), grid on;

% Tempo eliminazione farmaco 
tau = 1/min(abs(eig(A)));
te = 5*tau;

% Concentrazione minima e massima
c_min = min(ct);
c_max = max(ct);

% AUC 
auc2 = trapz(T,ct);

% Frazione
fa = k21/(k21+k01);

% Variazione k01 
k_var = linspace(k01/100,k01*100,10);
for i=1:10
    figure(5);
    % Matrici del sistema compartimentale lineare
    A_var = [-(k_var(i)+k21),0;k21,-k02]; 
    B = [d;0];
    C = [0,1/V2];
    D = 0;

    % Creazione del modello compartimentale lineare
    sys = ss(A_var,B,C,D);

    % Grafici della quantita e della concentrazione
    t = [0:0.01:10];
    [ct,T,qt] = impulse(sys,t);
    subplot(1,2,1), plot(T,qt), title("Quantita nel compartimento"), 
    xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), legend("Compartimento 1","Compartimento 2"), grid on;
    hold on;
    subplot(1,2,2), plot(T,ct), title("Concentrazione nel compartimento 2"),
    xlabel("Tempo (ore)"), ylabel("Concentrazione (mg/L)"), grid on;
    hold on;
end

% Calcolo della biodisponibilita
t_auc2 = [0.01, 0.23, 0.48, 1.55, 2.05, 3.01, 5.65];
c_auc2 = [2.15, 30.1323, 36.6604, 15.0529, 8.4495, 2.6962, 0.1136];

bio_disp = trapz(t_auc2,c_auc2) / auc1_nc;

% Calcolo costante assorbimento apparente (approccio non comp)
aumc2 = trapz(t_auc2,t_auc2.*c_auc2);
mrt2 = aumc2 / auc2;

mat = mrt2 - mrt;


%% Esercizio n*3 - Modello bicompartimentale (C-peptide)

% Dati modello 
k01 = 6.54e-2;
k12 = 5.68e-2;
k21 = 7.16e-2;
V1 = 3.29;
d = 49650;

% Matrici del sistema compartimentale lineare
A = [-(k01+k21),k12;k21,-k12]; 
B = [d;0];
C = [1/V1,0];
D = 0;

% Creazione del modello compartimentale lineare
sys = ss(A,B,C,D);

% Grafici della quantita e della concentrazione
figure(6);
t = [0:1:200];
[ct,T,qt] = impulse(sys,t);
subplot(1,2,1), plot(T,qt), title("Quantita nei compartimenti"), 
xlabel("Tempo (min)"), ylabel("Quantita (pmol)"), legend("Compartimento 1","Compartimento 2"), grid on;
subplot(1,2,2), plot(T,ct), title("Concentrazione nel compartimento 1"),
xlabel("Tempo (min)"), ylabel("Concentrazione (pmol/L)"), grid on;

% Grafico della concentrazione in scale logaritmica
figure(7);
semilogy(T,ct), title("Concentrazione nel compartimento 1 (scala log)"),
xlabel("Tempo"), ylabel("Concentrazione"), grid on;

% Tempo eliminazione C-peptide
tau_cp = 1/min(abs(eig(A)));
te_cp = 5*tau_cp; 

% AUC 
auc_cp = trapz(ct);

% CLtot 
CLtot = d/auc_cp;

% Variazione parametro
k_var = linspace(k21/100,k21*100,10);

figure(8);

for i=1:10
    % Matrici del sistema compartimentale lineare
    A_var = [-(k01+k_var(i)),k12;k_var(i),-k12];
    B = [d;0];
    C = [1/V1,0];
    D = 0;

    % Creazione del modello compartimentale lineare
    sys = ss(A_var,B,C,D);

    % Grafici della quantita e della concentrazione
    t = [0:1:200];
    [ct,T,qt] = impulse(sys,t);
    subplot(2,1,1), plot(T,qt), title("Quantita nei compartimenti"), 
    xlabel("Tempo (min)"), ylabel("Quantita (pmol)"), legend("Compartimento 1","Compartimento 2"), grid on;
    hold on;
    subplot(2,1,2), plot(T,ct), title("Concentrazione nel compartimento 1"),
    xlabel("Tempo (min)"), ylabel("Concentrazione (pmol/L)"), grid on;
    hold on;
end


%% Modello a tre compartimenti 

% Dati modello 
V1 = 5;
k21 = 2.22;
k12 = 0.859;
k31 = 0.031;
k13 = 0.008;
k01 = 1.2;
d = 500; 

% Matrici del sistema compartimentale lineare
A = [-(k01+k21+k31),k12,k13;k21,-k12,0;k31,0,-k13]; 
B = [d;0;0];
C = [1/V1,0,0];
D = 0;

% Creazione del modello compartimentale lineare
sys = ss(A,B,C,D);

% Grafici della quantita e della concentrazione
figure(9);
t = [0:1:200];
[ct,T,qt] = impulse(sys,t);
subplot(1,2,1), plot(T,qt), title("Quantita nei compartimenti"), 
xlabel("Tempo (min)"), ylabel("Quantita (pmol)"), legend("Compartimento 1","Compartimento 2","Compartimento 3"), grid on;
subplot(1,2,2), plot(T,ct), title("Concentrazione nel compartimento 1"),
xlabel("Tempo (min)"), ylabel("Concentrazione (pmol/L)"), grid on;

% Grafico della concentrazione in scale logaritmica
figure(10);
semilogy(T,ct), title("Concentrazione centrale (scala log)"),
xlabel("Tempo"), ylabel("Concentrazione");

% AUC 
auc_trc = trapz(ct);

% Volume di distribuzione
E = eig(A);
vd = d/(auc_trc*abs(min(E)));

% Pole zero map - Osservazione: quasi cancellazione zero-polo
figure(11);
pzplot(sys);

%% Dosi ripetute modello monocompartimentale

% Dati del modello 
ke = 1.2;
d = 500;
V1 = 5;

% Matrici del sistema compartimentale lineare
A = -ke;
B = d;
C = 1/V1;
D = 0;

% Creazione del modello compartimentale lineare
sys = ss(A,B,C,D);

% Dosi ripetute
deltaImpulsi = 6;
ampiezzaImpulsi = 1;
numeroImpulsi = 12;
t = linspace(0,deltaImpulsi,72);
xprec = 0;

figure(12);  

Y = [];
X = [];
T = [];

for i = 1:numeroImpulsi
    [yimp, timp, ximp] = impulse(sys*ampiezzaImpulsi,t);
    [yini,tini,xini] = initial(sys,xprec,t);
    
    Y = [Y;yimp+yini];
    X = [X;ximp+xini];
    T = [T,t+deltaImpulsi*(i-1)];
    
    xfinale=ximp+xini;
    xprec=xfinale(end,:);
    
end

plot(T,Y);
xlabel("Tempo (ore)"), ylabel("Quantita (mg)"), grid on;


figure(13);

% Delta = 3
[T1,Y1] = dosiRipetute(sys,3,1);
subplot(4,1,1), plot(T1,Y1), axis tight, grid on;
% Delta = 1
[T2,Y2] = dosiRipetute(sys,1,1);
subplot(4,1,2), plot(T2,Y2), axis tight, grid on;
% Delta = 0.5
[T3,Y3] = dosiRipetute(sys,0.5,1);
subplot(4,1,3), plot(T3,Y3), axis tight, grid on;
% Delta = 0.2
[T4,Y4] = dosiRipetute(sys,0.2,1);
subplot(4,1,4), plot(T4,Y4), axis tight, grid on;

% I parametri sono costante di eliminazione e intervallo tra una
% somministrazione e l'altra. 
