function [T,Y] = dosiRipetute(sys,deltaImpulsi,ampiezzaImpulsi)
    t = linspace(0,deltaImpulsi,72);
    numeroImpulsi = 72/deltaImpulsi;
    xprec = 0;  

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
end

