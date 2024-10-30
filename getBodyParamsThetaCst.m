function [mRb,mTb,theta] = getBodyParamsThetaCst(Lbody,dbody,eta,theta0)

    Nb = size(Lbody,1);
    p = Lbody./dbody;

    Zpara= 2*pi.*eta.*Lbody ./(log(p)-0.207+0.98./p-0.133./p.^2);
    Zperp = 4*pi.*eta.*Lbody ./(log(p)+0.839+0.185./p+0.233./p.^2);

    ZrPerp = pi/3*eta.*Lbody.^3./(log(p)-0.662+0.917./p-0.050./p.^2);
    ZrPara = (3.841*pi.*eta.*Lbody.*(dbody).^2)/4 .*(1+1.119 *10^(-4)+0.6884 ./p -0.2019./p.^2); % fit Tirado 1980

    % Effective friction coefficients
    
    theta = theta0;
    
    mRb = ZrPerp.*sin(theta).^2 + ZrPara.*cos(theta).^2;
    
    mTb = (Zpara + (Zperp-Zpara).* sin(theta).^2);

end