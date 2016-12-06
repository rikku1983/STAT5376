function Fstr = regfngrad(xfine, y0fine, dhdc, yregfd, Wfd, ...
                          Kmat, periodic, crit, nvar)
  
  nfine   = length(xfine);
  cvec    = getcoef(Wfd);
  ncvec   = length(cvec);
  onecoef = ones(1,ncvec);
  
  if periodic
     dhdc(:,1) = 1;
  else
     dhdc(:,1) = 0;  
  end
  yregmat  = squeeze(eval_fd(xfine, yregfd));
  Dyregmat = squeeze(eval_fd(xfine, yregfd, 1));
  
  %  loop through variables computing function and gradient values
  
  Fval = 0;
  gvec = zeros(ncvec,1);
  for ivar = 1:nvar
    y0ivar  =   y0fine(:,ivar);
    ywrthi  =  yregmat(:,ivar);
    Dywrthi = Dyregmat(:,ivar);
    aa = mean(y0ivar.^2);
    cc = mean(ywrthi.^2);
    bb = mean(y0ivar.*ywrthi);
    ff = aa - cc;
    dd = sqrt(ff^2 + 4*bb^2);
    Dywrtc  = (Dywrthi * onecoef).*dhdc;
    if crit == 1
        %  least squares criterion
      res  = y0ivar - ywrthi;
      Fval = Fval + aa - 2*bb + cc;
      gvec = gvec - 2.*Dywrtc'*res./nfine;
    elseif crit == 2
         % minimum eigenvalue criterion
      Fval = Fval + aa + cc - dd;
      Dbb  =    Dywrtc'*y0ivar./nfine;
      Dcc  = 2.*Dywrtc'*ywrthi./nfine;
      Ddd  = (4.*bb.*Dbb - ff.*Dcc)./dd;
      gvec = gvec + Dcc - Ddd;
    elseif crit == 3
        %  least squares with root Dh weighting criterion
      res  = y0ivar - ywrthi.*sqrt(dhdt);
      Fval = Fval + mean(res.^2);
      gvec = gvec - 2.*Dywrtc'*res./nfine;
    else
        error('Inadmissible value for CRIT.');
    end
  end
  if ~isempty(Kmat)
     ind   = 2:ncvec;
     ctemp = cvec(ind,1);
     Kctmp = Kmat*ctemp;
     Fval  = Fval + ctemp'*Kctmp;
     gvec(ind) = gvec(ind) + 2.*Kctmp;
  end
  
%  set up F structure containing function value and gradient

  Fstr.f    = Fval;
  Fstr.grad = gvec;
  %  do not modify initial coefficient for B-spline and Fourier bases
  if ~periodic,  Fstr.grad(1) = 0;  end
  Fstr.norm = sqrt(sum(Fstr.grad.^2));
