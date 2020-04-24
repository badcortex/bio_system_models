function diff = funzionale(p,t,z,m)
    switch m 
        case 1
            diff = z - p(1)*exp(-p(2)*t);
        case 2    
            diff = z - (p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t));
        case 3
            diff = z - (p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t) + p(5)*exp(-p(6)*t));
    end
end

