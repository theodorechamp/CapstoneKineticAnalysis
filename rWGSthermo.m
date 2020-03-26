function output = rWGSthermo()
%RWGSTHERMO calculates the thermodynamic equilibrium at various
%temperatures and plots the results.
    output = 0; %#ok<NASGU> %if code fails, return 0
    
%all temperatures are in Kelvin
    Trange = linspace(373,1873,200);
    l = length(Trange);

    R = 8.314; %J/mol-K or L-kPa/mol-K
    Tref = 298; %K

    dGref1 = -28519.6; %J/mol
    dGref2 = -141941.979; %J/mol
    Keqref1 = exp(-dGref1/(R*Tref));
    Keqref2 = exp(-dGref2/(R*Tref));

    dHrxn1 = -41000; %J/mol, assume constant with temp
    dHrxn2 = -206000; %J/mol, assum constant with temp

%calculate Keq at each temperature
    Keq = zeros(l,2);
    n0 = [0,0,1,0,3];% n = [nH2, nCO, nCO2, nH2O, nCH4]
    V = 20; %L
    
    for i = 1:l
        Keq1 = Keqref1 * exp(-dHrxn1/R*(1/Trange(i)-1/Tref));
        Keq2 = Keqref2 * exp(-dHrxn2/R*(1/Trange(i)-1/Tref));
%         P = sum(n0) * R * Trange(i)/V;
        P = 100;
        Keq2 = Keq2*P^2; % account for pressure change from reaction
        Keq(i,:) = [Keq1 Keq2];
    end

%calculate compositions at each temperature
    comp = zeros(l,5); %comp = [PH2, PCO, PCO2, PH2O, PCH4]
    x0 = [0,0];
    
     A = [ -1,  3;
            1,  1;
           -1,  0;
            1, -1;
            0, -1 ];
     b = [n0(1);n0(2);n0(3);n0(4);n0(5)];
     options = optimoptions('fgoalattain','MaxFunctionEvaluations',1000, 'display', 'notify');

     for i = 1:l
        x = fgoalattain(@(x)rWGSKeqs(x,n0,Keq(i,:)),x0,[0,0],[0.001,0.001],A,b,[],[],[],[],[],options);
        
        comp(i,1)  = n0(1) +   x(1) - 3*x(2);
        comp(i,2)  = n0(2) -   x(1) -   x(2);
        comp(i,3)  = n0(3) +   x(1) -      0;
        comp(i,4)  = n0(4) -   x(1) +   x(2);
        comp(i,5)  = n0(5) +      0 +   x(2);
        comp(i,:) = comp(i,:)/sum(comp(i,:));
        %comp(i,6:7)  = x(1,:);
    end

%plot the results
    close all
    hold on
    figure(1)
    plot(Trange,comp(:,1:5))
    legend('H2','CO','CO2','H2O','CH4')
    
    output = [Trange' Keq comp];
end

