%IVGTT: IntraVenous Glucouse Tolerance Test
% OGTT: Oral Glucose Tolerance Test

% Sono test in cui uno stimolo di una certa natura viene usato per
% controllare la risposta del paziente e capire se esula da un'andamento
% "normale"

% In questa esercitazione, l'Up&Down è uno stimolo(glucosio) ad andamento 
% variabile nel tempo e si raccolgono campioni, in plasma, che siano 
% indicativi della secrezione di insulina.
% Abbiamo perciò i dati di uscita (valore insulina nel compartimento
% accessibile) e vogliamo ricostruire l'ingresso.
% In realtà, come gia' detto in precedenza, si misura l'andamento di
% C-Peptide, la cui dinamica è lineare nel range di valori di studio e
% viene prodotto in rapporto equimolare con l'insulina.

clear all
close all
clc

%ESERCITAZIONE 7
%abbiamo salvato in un file datiesrc7 la stima dei beta per ogni
%paramentro globale

%inoltre, utilizzando i dati del soggetto forniti dal testo abbiamo
%calcolato il valore dei paramentri con l'esercitazione 6 lanciata per ogni
%paramentro globale


%carico i dati del file datieserc7 e li salvo in nuove variabili
datiEs7;    %tempo, conc e basale
stimaEs7; %nostri beta

minuti=DATA(:,1); %tempi
concentrazione=DATA(:,2); %concentrazioni
basale=yb;


% vol=3.8357;
% temp_emi_short=4.9963;
% temp_emi_long=34.5287;
% fraction=0.76502;
%--------------------------------------------------------------------------
%Applico la deconvoluzione discreta cioè
%ipotizzo di ricostruire u come una funzione costante a tratti per
%risolvere il fatto che la deconvoluzione è un problema mal posto

%RAW DECONVOLUTION : NO RUMORE

%data l'uscita g(solitamente somma di esponenziali) e i dati di
%concentrazione vorrei ricostruire l'ingresso.
%y=G*u
%in assenza di rumore dovrei fare G^-1*z e otterrei esattamente u
%Per poter invertire G--> G QUADRATA, cioè significa che la griglia
%virtuale (dell'ingresso) coincide con quella di campionamento (dell'uscita)data dal file letto
%Ipotizzo di ricostruire u come una funzione costante a tratti per
%risolvere il fatto che la deconvoluzione è un problema mal posto

%costruisco la tipica uscita g(t)=A*exp(-alpha)+B*exp(-beta)

 Alpha = log(2)/temp_emi_short;
 %il tempo di emivita di distribuzione è log(2)/alfa
 Beta = log(2)/temp_emi_long;
 %il tempo di emivita è log(2)/beta
 A = fraction/vol;
 B = (1 - fraction)/vol;
 
 n=length(concentrazione);
 minuti0=[0 minuti']; % perche se no la trapz non considera la prima parte di area (tra il minuto 0 e il minuto 10)
                      % le misurazioni nel set di dati partono dal minuto
                      % 10
 
 for i=1:length(minuti)

   intervallo=  [minuti0(i) minuti0(i+1)];     % intervallo sotto il quale la trapz calcola l'area
   %per per precauzione, siccome potrebbe capitare il caso di un campionamento non uniforme, calcolo le aree per
   %ogni intervallo temporale tra i campioni, sostanzialmente ho un
   %intervallo ogni dieci minuti tot=24
   g(i)=trapz(intervallo,A*exp(-Alpha*intervallo)+B*exp(-Beta*intervallo));
 
 end
 
 G=zeros(n); %nxn
 
 %COSTRUZIONE G
 
 for i=1:n %colonne
    for j=i:n %righe
       G(j,i)=g(j-i+1);  %compilo colonna per colonna costruendo una matrice triangolare 
                         %inferiore che 
                         %rappresenta lo scivolare nel tempo (verso dx) della funzione 
                         %che esce dalla finestra temporale
        end
    end

z=concentrazione-basale;     %variazione di concentrazione dal basale z=y-yb

u_hat=G^-1*z;        %vettore degli ingressi stimati con la raw deconvolution u=G^-1*z

z_riconvoluto=G*u_hat;     %riconvoluzione partendo dagli ingressi stimati

residui = z - z_riconvoluto;    %differenza tra le concentrazioni effettive e quelle stimate 
                                %con la riconvoluzione
                                %nel caso della raw convolution si hanno
                                %per definizione i residui=0
 

%GRAFICI
%u_hat è costante a tratti, quindi devo plottarla in modo che il suo valore
%rimanga costante nell'intervallo di tempo in cui considero il singolo
%valore di u_hat stesso: uso l'istruzione stairs

figure(1)

subplot(3, 1, 1)
U_hat=[u_hat', u_hat(end)]; %duplico ultimo valore di u_hat per costruire 
                            %bene l'ultima parte del grafico
stairs(minuti0, U_hat, 'r');
xlabel('tempi'); ylabel('u hat')
hold on
plot(minuti, u_hat, '*')
title('Raw-deconvolution')


subplot(3,1,2)
plot(minuti, z_riconvoluto, minuti, z, 'og')
xlabel('tempi'); ylabel('concentrazioni')
title('Riconvoluzione e dati originari')


subplot(3,1,3)
plot(minuti, residui, 'og',minuti,zeros(length(minuti),1),'r')
xlabel('tempi'); ylabel('residui')
ylim([-0.002 0.002]) %scala ragionevole per osservare l'andamento nullo dei residui in questo caso
title('Residui')

%PROBLEMI RAW DECONVOLUTION
%malcondizionamento: una minima oscillazione della y viene amplificata
%tantissimo sull'ingrsso deconvoluto
%griglia virtuale= griglia di campionamente (vorrei virtuale più fitta e 
%campionamento meno per evitare più misure)
%assunzione che i dati fossero privi di rumore

%--------------------------------------------------------------------------
%REGOLARIZZAZIONE : CADE IPOTESI DI ASSENZA DI RUMORE
%y>>numero incognite
%infittisco la griglia virtuale ma il problema risulta ancora risolvibile
%perchè gli ingressi NON sono INDIPENDENTI tra loro: sto introducendo l'hp
%che l'ingresso deconvoluto abbia una dinamica poco rapida 

griglia_virtuale = [min(minuti0): 1 : max(minuti0)]; %va da 0 a 240 con passo 1--> 
                                                     %costruisco griglia temporale 
                                                     %da 0 a 240
m=length(griglia_virtuale); %241


 for k=1:m-1
   intervallo_reg= [griglia_virtuale(k) griglia_virtuale(k+1)];   %la prima volta è [0 1] 
   
   g_reg(k)=trapz(intervallo_reg,A*exp(-Alpha*intervallo_reg)+B*exp(-Beta*intervallo_reg));
 
 end
 
 G_reg=zeros(m-1); %240x240
 
 %COSTRUZIONE G
 
 for k=1:(m-1)
    for l=k:(m-1)
      G_reg(l,k)=g_reg(l-k+1);  %compilo colonna per colonna (triangolare inferiore)
    end
 end

    
 
%CREO LE MATRICI CHE MI SERVONO NEI CRITERI PER LA RICERCA DEL GAMMA OTTIMO
%Trovo gli indici degli istanti di campionamento corrispondenti agli istanti
%della griglia virtuale e considero solo queste righe.


for i=1:length(minuti)
    indici(i)=find(minuti(i)==griglia_virtuale); %trovo gli indici corrispondenti 
                                                 %alla vecchia griglia (10 20 30 40...)
end

G_reg_new=G_reg(indici-1, :); %prendo le 24 righe corrispondenti alla vecchia griglia e le metto in una matrice G_reg_new
                              %la matrice G_REG in questo caso non deve
                              %essere quadrata, infatti sarà del tipo 24x24

%costruiamo la matrice penalità considerando la derivata prima (K=1)
%P è la somma di 2 matrici:
%il primo termine è una matrice nulla tranne la
%sottodiagonale (-1 per tutti gli elementi)
%il secondo termine è una matrice identità.

P=diag(-ones(1,length(g_reg)-1),-1) + diag(ones(1,length(g_reg))); %240x240

%costruiamo la matrice sigmaV (VARIANZA) considerando un CV costante e pari a 0.04
%(preso dalle altre esercitazioni)

CV=0.04;
mat=diag(concentrazione.^2);
sigmaV=CV*CV*mat; %sigmaV=CV^2*matrice diagonale con conc al quadrato


%% RICERCA DEL GAMMA OTTIMO.

%Esistono 3 criteri di ricerca del gamma ottimo:


% 1)CRITERIO DI DISCREPANZA 

%Calcolo il gamma per il quale le caratteristiche statistiche dei residui
%sono uguali a quelli dell'errore di misura 

%Dal punto di vista pratica mi interessa che la somma al quadrato dei
%residui sia uguale alla traccia della matrice sigma_v


gamma=0.01;  %si parte da un gamma molto piccolo (riproduzione fedele dei dati) 
             %per poi incremetarlo progressivamente
             %più gamma è grande, più è forte la regolarizzazione

             
%calcolo la traccia di sigma_v
traccia= trace(sigmaV);

%Stimo l'ingresso a partire dall'uscita
u_reg=(inv(G_reg_new'*inv(sigmaV)*G_reg_new+gamma*P'*P))*G_reg_new'*inv(sigmaV)*z; 

%u_hat=((G'*sigmaV^-1*G)+Gamma*P'*P)^-1*(G'*sigmaV^-1*uscite)
%Stima dell'uscita dati gli ingressi stimati  
y_rec=G_reg_new*u_reg;

residui=z - y_rec;
RSS=norm(residui)^2; %Somma dei quadrati dei residui 

%modifico il gamma finchè la somma dei quadrati dei residui=somma varianze
%degli errori di misura--> somma|e'e|=E[v'v]=Somma(sigmaV)=traccia(sigma V).

while (abs(RSS-traccia)>0.04) 
    gamma=gamma*2;
    u_reg=(inv(G_reg_new'*inv(sigmaV)*G_reg_new+gamma*P'*P))*G_reg_new'*inv(sigmaV)*z;
    y_rec=G_reg_new*u_reg;
    residui=z - y_rec;
    RSS=norm(residui)^2;
    
end

gamma_Twomey=gamma;

%SVANTAGGI: mancanza di solide basi teoriche, rischio di oversmoothing


% 2)CRITERIO GVC:Generalized Cross Validation.
%Gamma calcolato come argmin(funzionale di costo GVC) 

%Dato WRSS=somma pesata dei residui al quadrato
%GCV=WRSS/(TRACCIA(In-psi))^2
%Dove psi dipende da Gamma,G,SigmaV
%calcolato nella funzione metodo_gvc

%Facciamo partire Gamma da un valore alto(10^3)

gamma_GVC=fminsearch(@metodo_gvc,1000,[],z,sigmaV,G_reg_new,P)

[GVC u_GVC y_GVC residui_GVC]=metodo_gvc(gamma_GVC,z,sigmaV,G_reg_new,P);

% 3)CRITERIO ML

%Scelgo come valore iniziale di gamma quello trovato con il primo
%criterio(discrepanza)

gamma=gamma_Twomey;

u_ML=(inv(G_reg_new'*inv(sigmaV)*G_reg_new+gamma*P'*P))*G_reg_new'*inv(sigmaV)*z;
y_ML=G_reg_new*u_ML;
residui=z - y_ML;

psi_ML= G_reg_new * inv(G_reg_new'*(sigmaV)^-1*G_reg_new + gamma*P'*P) * G_reg_new'*(sigmaV)^-1; %matrice di influeza
q=trace(psi_ML);  %gradi di libertà
gamma_ML= (residui'*(sigmaV^-1)*residui*q)/(n-(q*u_ML'*P'*P*u_ML));

valore1=gamma;
valore2=gamma_ML;

while (abs(valore2-valore1) > 0.04) %Criterio di taratura (tolleranza)
    
    valore1=valore2;
    u_ML=(inv(G_reg_new'*inv(sigmaV)*G_reg_new+gamma_ML*P'*P))*G_reg_new'*inv(sigmaV)*z;
    y_ML=G_reg_new*u_ML;
    residui=z - y_ML;
    psi_ML= G_reg_new * inv(G_reg_new'*(sigmaV)^-1*G_reg_new + gamma_ML*P'*P) * G_reg_new'*(sigmaV)^-1;
    q=trace(psi_ML);
    gamma_ML= (residui'*sigmaV^-1*residui*q)/(n-(q*u_ML'*P'*P*u_ML));
   
    valore2=gamma_ML;
   
end


%GRAFICI  
figure(2)

% si plottano i tre criteri su 3 colonne

%%PER OGNI CRITERIO : INGRESSO STIMATO,USCITA RICONVOLUTA,RESIDUI

subplot(3,3,1)
tempi_u_reg = griglia_virtuale(1:end-1);  % creo il vettore dei tempi per u_reg considerando i valori della griglia virtuale (escluso il valore finale)
plot(tempi_u_reg,u_reg)
hold on
stairs(minuti0,U_hat, 'r')
xlabel('tempo');ylabel('u_reg')
title(['Regolarizzazione: gamma Twomey= ' num2str(gamma)])

subplot(3,3,2)
plot(minuti, y_rec, minuti, z, 'og')
xlabel('tempo'); ylabel('concetrazioni predette')
title('Riconvoluzione e dati')

subplot(3,3,3)
plot(minuti, residui, 'og', minuti, residui, ':b')
xlabel('tempo'); ylabel('y - yhat')
title('Residui')

subplot(3, 3, 4)
plot(tempi_u_reg,u_GVC)
hold on
stairs(minuti0, U_hat, 'r')
xlabel('tempo'); ylabel('u_reg')
title(['Regolarizzazione: gamma GVC= ' num2str(gamma_GVC)])
axis([0 250 -0.1 1])

subplot(3,3,5)
plot(minuti, y_GVC, minuti, z, 'og')
xlabel('tempo'); ylabel('concetrazioni predette')
title('Riconvoluzione e dati')

subplot(3,3,6)
plot(minuti, residui_GVC, 'og', minuti, residui_GVC, ':b')
xlabel('tempo'); ylabel('y - yhat')
title('Residui')

subplot(3, 3, 7)
plot(tempi_u_reg,u_ML)
hold on
stairs(minuti0, U_hat, 'r')
xlabel('tempo'); ylabel('u_ML')
title(['Regolarizzazione: gamma ML= ' num2str(gamma_ML)])
axis([0 250 -0.1 1])

subplot(3,3,8)
plot(minuti, y_ML, minuti, z, 'og')
xlabel('tempo'); ylabel('concentrazioni predette')
title('Riconvoluzione e dati')

subplot(3,3,9)
plot(minuti, residui, 'og', minuti, residui, ':b')
xlabel('tempo'); ylabel('y - yhat')
title('Residui')


% I tre metodi portano a valori di gamma molto diversi:

%-criterio di Twomey: alto (ordine di 10^3)
%-criterio GVC: intermedio (ordine di 10^2)
%-criterio ML: molto basso (ordine di 10^-8)

%Siccome il gamma è¨ il parametro di regolarità , più¹ gamma è¨ alto più u è regolare
%e i gradi di libertà diminuiscono. 





