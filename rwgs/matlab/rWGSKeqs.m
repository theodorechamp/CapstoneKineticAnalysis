function F = rWGSKeqs(x,n,Keq)
%function F = rWGSKeqs(x,n, Keq)
%   n = [nH2, nCO, nCO2, nH2O, nCH4]
    
    nH2  = n(1) +   x(1) - 3*x(2);
    nCO  = n(2) -   x(1) -   x(2);
    nCO2 = n(3) +   x(1) -      0;
    nH2O = n(4) -   x(1) +   x(2);
    nCH4 = n(5) +      0 +   x(2);
    ntot = nH2+nCO+nCO2+nH2O+nCH4;

    F(1) = Keq(1)*(nCO*nH2O) - (nCO2*nH2);
    F(2) = Keq(2)*(nCO*nH2^3) - (nCH4*nH2O)*(ntot)^2;
end
   
