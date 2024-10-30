function [mRb,mTb,theta] = getBodyParamsComputeTheta(Lbody,dbody,eta,Nf,mTT,mTR,mRR)

    Nb = size(Lbody,1);
    p = Lbody./dbody;

    Zpara= 2*pi.*eta.*Lbody ./(log(p)-0.207+0.98./p-0.133./p.^2);
    Zperp = 4*pi.*eta.*Lbody ./(log(p)+0.839+0.185./p+0.233./p.^2);

    ZrPerp = pi/3*eta.*Lbody.^3./(log(p)-0.662+0.917./p-0.050./p.^2);
    ZrPara = (3.841*pi.*eta.*Lbody.*(dbody).^2)/4 .*(1+1.119 *10^(-4)+0.6884 ./p -0.2019./p.^2); % fit Tirado 1980

    % Effective friction coefficients
    
    mSin = zeros(Nb,length(Nf));
    for n=1:length(Nf)
        alpha = 2*pi*rand(Nb,Nf(n));

        Ssin = sum(sin(alpha),2);
        Scos = sum(cos(alpha),2);

        a0 = -atan(Scos./Ssin);
        
        mSin(:,n) = mean(sin(alpha-a0) ,2);
        
        
    end
    muTaprx = (Zpara+Zperp)/2;
    mRed = mTR./(mRR - mTR.^2./(mTT+muTaprx)) .* muTaprx./(mTT+muTaprx);
    theta = abs( 0.5*asin(2* mSin.*ZrPara./(ZrPerp-ZrPara).* dbody/2 .* mRed)  );
    
    
    mRb = ZrPerp.*sin(theta).^2 + ZrPara.*cos(theta).^2;
    
    mTb = (Zpara + (Zperp-Zpara).* sin(theta).^2);

end