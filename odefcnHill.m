function dqdt = odefcnHill(t,q,e,k01,k02,Vmax,km)
    dqdt = zeros(2,1);
    dqdt(1) = -(k01+Vmax*(q(1).^(e-1))/(km^e+(q(1).^e)))*q(1);
    dqdt(2) = (Vmax*(q(1).^(e-1))/(km^e+(q(1).^e)))*q(1)-k02*q(2);
end

