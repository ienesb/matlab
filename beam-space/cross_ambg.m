function C = cross_ambg(taus, nus, T)   
    % taus is row vector, nus is column vector, T is scalar
    C = zeros(length(taus), length(nus));
    
    C(taus<0, nus==0) = T + repmat(taus(taus<0),1,sum(nus==0));
    C(taus>=0, nus==0) = T - repmat(taus(taus>=0),1,sum(nus==0));
    
    C(taus<0, nus~=0) = 1j*(exp(-1j*2*pi*(T+taus(taus<0))*nus(nus~=0))-1)./(2*pi*nus(nus~=0));
    C(taus>=0, nus~=0) = 1j*(exp(-1j*2*pi*nus(nus~=0)*T)-exp(-1j*2*pi*taus(taus>=0)*nus(nus~=0)))./(2*pi*nus(nus~=0));
    
    C(taus<-T, :) = 0; 
    C(taus>T, :) = 0; 

    C = C ./ T;
end
