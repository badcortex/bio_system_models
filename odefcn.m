function dqdt = odefcn(t,q,k01,k02,Vmax,km)
    dqdt = zeros(2,1);
    dqdt(1) = -(k01+Vmax/(km+q(1)))*q(1);
    dqdt(2) = (Vmax/(km+q(1)))*q(1)-k02*q(2);
end

