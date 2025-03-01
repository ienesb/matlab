function H = getH(n, np, m, mp, chanParams, T) 
    df = 1/T;

    H = 0;
    for p = 1:3
        rhop = chanParams.pathGains(p);
        taup = chanParams.pathDelayTimes(p);
        nup = chanParams.pathDopplerFreqs(p);
        C = cross_ambg((n-np)*T-taup, (m-mp)*df-nup, T);
    
        H = H + rhop*exp(1j*2*pi*taup*nup)*C*exp(1j*2*pi*(np*T*nup-mp*df*taup));
    end

end