function [Data,Fct] = generateModelFlagellarPropulsion(model,Lflag,Nflag,lmbd0,R0,r0,Lb,db,varargin)
    gkCst = 0.5;
    eta = 1; % mPa.s
    
    % parse data
    pars = inputParser;
    addOptional(pars,'dlmbd',0);
    addOptional(pars,'dR',0);
    addOptional(pars,'dr0',0);
    addOptional(pars,'dL',0);
    addOptional(pars,'dd',0);
    addParameter(pars,'gkcst',gkCst,@isnumeric);
    addParameter(pars,'Nstat',1,@isnumeric);
    addParameter(pars,'wm',220,@isnumeric);
    addParameter(pars,'eta',eta,@isnumeric);
    addParameter(pars,'theta0','none',@(x) isnumeric(x)||strcmp(x,'none'));
    parse(pars,varargin{:});
    dlmbd = pars.Results.dlmbd;
    dR = pars.Results.dR;
    dr = pars.Results.dr0;
    dL = pars.Results.dL;
    dd = pars.Results.dd;
    gkCst = pars.Results.gkcst;
    n = pars.Results.Nstat;
    wM = 2*pi*pars.Results.wm;
    eta = pars.Results.eta;
    theta0 = pars.Results.theta0;
    
    % set functions
    switch model
        case 'Lighthill'
            ratio = @(lmbd,psi,r,N) lmbd./cos(psi)./(r.*sqrt(N));
            kn = @(lmbd,psi,r,N) 4*pi*eta./(log(0.18*ratio(lmbd,psi,r,N))+0.5); 
            gk = @(lmbd,psi,r,N) 0.5*(log(0.18*ratio(lmbd,psi,r,N))+0.5)./(log(0.18*ratio(lmbd,psi,r,N)));
        case 'GrayHancock'
            ratio = @(lmbd,psi,r,N) lmbd./(r.*sqrt(N));
            kn = @(lmbd,psi,r,N) 4*pi*eta./(log(2*lmbd./(r.*sqrt(N)))+0.5); 
            gk = @(lmbd,psi,r,N) 0.5*(log(2*lmbd./(r.*sqrt(N)))+0.5)./(log(2*lmbd./(r.*sqrt(N)))-0.5);
        case 'Lighthill-Constantgk'
            ratio = @(lmbd,psi,r,N) lmbd./cos(psi)./(r.*sqrt(N));
            kn = @(lmbd,psi,r,N) 4*pi*eta./(log(0.18*ratio(lmbd,psi,r,N))+0.5); 
            gk = @(lmbd,psi,r,N) gkCst;
        case 'GrayHancock-Constantgk'
            ratio = @(lmbd,psi,r,N) lmbd./(r.*sqrt(N));
            kn = @(lmbd,psi,r,N) 4*pi*eta./(log(2*lmbd./(r.*sqrt(N)))+0.5);  
            gk = @(lmbd,psi,r,N) gkCst;
        case 'RC'
            ratio = @(lmbd,psi,r,N) lmbd./(r.*sqrt(N));
            kn = @(lmbd,psi,r,N) 4*pi*eta./(log(0.25*lmbd./(r.*sqrt(N)))+0.5); 
            gk = @(lmbd,psi,r,N) 0.7;
        otherwise
            ex = MException('Unsupported model');
            throw(ex)
    end
    
    mTT = @(l,lmbd,R,psi,r,N) kn(lmbd,psi,r,N) .* l.*sin(psi).*(tan(psi)+gk(lmbd,psi,r,N)./tan(psi));  
    mTR = @(l,lmbd,R,psi,r,N) kn(lmbd,psi,r,N) .* l.*R.*sin(psi).*(1-gk(lmbd,psi,r,N));
    mRR = @(l,lmbd,R,psi,r,N) kn(lmbd,psi,r,N) .* l.*R.^2.*sin(psi).*(gk(lmbd,psi,r,N).*tan(psi)+1./tan(psi));
    
    Fct.ratio = ratio;
    Fct.kn = kn;
    Fct.gk = gk;
    Fct.mTT = mTT;
    Fct.mTR = mTR;
    Fct.mRR = mRR;
    
    % Prepare parameters
    % flagellum
    lmbd = (lmbd0 + dlmbd * (rand(n,1)-0.5) );
    R = (R0 + dR * (rand(n,1)-0.5) );
    r = (r0 + dr * (rand(n,1)-0.5) );
    psi = (atan((2*pi*R)./lmbd) );
    
    Lf = reshape(Lflag, 1, length(Lflag(:)));
    Nf = reshape(Nflag, 1, length(Nflag(:)));
    
    % body
    Lbody = Lb + dL * (rand(n,1)-0.5);
    dbody = db + dd * (rand(n,1)-0.5);
    
    % Compute coefficients
    mTTf = mTT(Lf,lmbd,R,psi,r,Nf);
    mTRf = mTR(Lf,lmbd,R,psi,r,Nf);
    mRRf = mRR(Lf,lmbd,R,psi,r,Nf);
    
    if ischar(theta0)
        [mRb,mTb,theta] = getBodyParamsComputeTheta(Lbody,dbody,eta,Nf,mTTf,mTRf,mRRf);
    else
        [mRb,mTb,theta] = getBodyParamsCorrectedv3(Lbody,dbody,eta,theta0);
    end
        
    
    mRRred = mRR(Lf,lmbd,R,psi,r,Nf)-mTR(Lf,lmbd,R,psi,r,Nf).^2./(mTT(Lf,lmbd,R,psi,r,Nf)+mTb);
    U_wf = 2*pi* mTR(Lf,lmbd,R,psi,r,Nf)./(mTT(Lf,lmbd,R,psi,r,Nf) + mTb);
    wb_wf = mRRred./mRb;
    U = mTR(Lf,lmbd,R,psi,r,Nf)./(mTT(Lf,lmbd,R,psi,r,Nf) + mTb) .* mRb./ (mRb + mRRred) .* wM;
    
    GammaM = mRb .* wb_wf ./(1+wb_wf) .* wM ./Nf;
    
    % Export data
    
    Data.mean.mRRred = reshape(mean(mRRred,1),size(Lflag));
    Data.mean.U_wf = reshape(mean(U_wf,1),size(Lflag));
    Data.mean.wb_wf = reshape(mean(wb_wf,1),size(Lflag));
    Data.mean.U = reshape(mean(U,1),size(Lflag));
    Data.mean.GammaM = reshape(mean(GammaM,1),size(Lflag));
    
    Data.std.mRRred = reshape(std(mRRred,[],1),size(Lflag));
    Data.std.U_wf = reshape(std(U_wf,[],1),size(Lflag));
    Data.std.wb_wf = reshape(std(wb_wf,[],1),size(Lflag));
    Data.std.U = reshape(std(U,[],1),size(Lflag));
    Data.std.GammaM = reshape(std(GammaM,[],1),size(Lflag));
    
    Data.mean.mRb = mean(mRb,1);
    Data.mean.mTb = mean(mTb,1);
    Data.mean.theta = mean(theta,1);
    
    Data.std.mRb = std(mRb,[],1);
    Data.std.mTb = std(mTb,[],1);
    Data.std.theta = std(theta,[],1);
    
    Data.stats = [max(Lbody),min(Lbody),max(dbody),min(dbody)];
    
    Data.model = model;

end