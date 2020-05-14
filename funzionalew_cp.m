function [residual,func_evals] = funzionalew_cp(p,t,z,w,d)
    % Matrici del sistema LTI in funzione dei parametri 
    A = [-(p(1)+p(3)),p(2);p(3),-p(2)];
    B = [d;0];
    C = [1/p(4),0];
    D = 0;
    
    % Creazione del sistema e risposta impulsiva
    sys = ss(A,B,C,D);
    ts = t(1):1:t(end);
    y = impulse(sys,ts);
    [~,ia,~] = intersect(ts,t);
    
    func_evals = y(ia);
    residual = w*(z-func_evals);
end

