function diff = funzionalew(p,t,z,m,w)
    switch m 
        case 1
            diff = w*(z - p(1)*exp(-p(2)*t));
        case 2    
            diff = w*(z - (p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t)));
        case 3
            diff = w*(z - (p(1)*exp(-p(2)*t) + p(3)*exp(-p(4)*t) + p(5)*exp(-p(6)*t)));
    end
end



