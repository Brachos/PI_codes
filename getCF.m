function Cf = getCF(length,M)
    data0 = load("cf0.csv");
    data07 = load("cf07.csv");
    data1 = load("cf1.csv");
    
    a = 303.1; %30000ft
    v = M*a; %vitesse
    u = 1.458*10^-5; %dynamic viscosity

    Re = v * length / u;
    
    cfM = [0,0.7,1];
    cf = [0,0,0];
    cf(1) = interp1(data0(:,1),data0(:,2),Re);
    cf(2) = interp1(data07(:,1),data07(:,2),Re);
    cf(3) = interp1(data1(:,1),data1(:,2),Re);
    
    Cf = interp1(cfM,cf,M);
end