function [gcv,u,y,res] = gvc_regularization(gamma,z,sigma_v,G,P)

    u = (inv(G'*inv(sigma_v)*G+gamma*P'*P))*G'*inv(sigma_v)*z;
    y = G*u;
    res = z - y;

    n = length(z); 

    % Calcolo della matrice di influenza
    H = G*inv(G'*inv(sigma_v)*G+gamma*P'*P)*G'*inv(sigma_v); 
    
    % Gradi di liberta equivalenti
    q = trace(H); 
    WRSS = (res'*inv(sigma_v)*res);
    gcv = WRSS/(n-q)^2; 
end