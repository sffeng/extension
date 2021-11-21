function Z=plot_turning
  rng(1985)
  close all
  n = 50; 
  d = 2;
  tv = [0 1 2]; % Time lags

  figure
  for s = 1:3
    subplot(2,3,s)
    t=tv(s);
    plothdxt(n,t,d)
    hold on
    plotunitcircle()
    axis square
  end
  subplot(2,3,1)
  ylabel('Covariance')


  %% now generate a single realization for all three time lags
  S = [];
  mesh = [];
  for s = 1:3
    t = tv(s);
    [X,Y,SS] = getcorrS(n,t,d);
    mm = [X(:) Y(:) repmat(t,n^2,1)];
    
    mesh = [mesh; mm];
    S = [S; SS(:)];
  end
  corr.name = 'gauss';
  corr.c0 = [1];
  corr.sigma=S;
  [F,KL] = randomfield(corr,mesh,'trunc',15);
  Z = reshape(F,n,n,3);

  for s = 4:6
    subplot(2,3,s)
    surf(X,Y,Z(:,:,s-3))
    hold on
    plotunitcircle()  
    view(2)
    shading interp
    if s==4
      ca = caxis;
    end
    caxis(ca)
    hold on
    plotunitcircle()
    axis square
  end

  subplot(2,3,4)
  ylabel('Realization')

  
  return
  
  function plotunitcircle()
    th = 0:.01:2*pi;
    one=ones(size(th));
    plot3(cos(th),sin(th),10*one,'linewidth',2,'color','k')
  
  function plothdxt(n,t,d)
    [X1,X2,Z] = getcorrS(n,t,d);
  surf(X1,X2,Z)
  shading interp
  caxis([0.3,.8])
  view(0,90)
  title(['t=' num2str(t)])
%  xlabel('x_1')
%  ylabel('x_2')
%  zlabel('h_d')
  
  function [X1,X2,S] = getcorrS(n,t,d)
  xx = linspace(-2,2,n);
  [X1,X2] = meshgrid(xx,xx);
  Z = zeros(n);
  for j = 1:n
    for k = 1:n
      x1 = X1(j,k);
      x2 = X2(j,k);
      X = norm([x1,x2],2);
      Z(j,k) = hdxt(X,t,d);
    end
  end
  S = Z;
      
  function val = hdxt(x,t,d)
    ff = @(u) hbar(u,t) .* (1-(u.^2)./(x.^2)).^((d-3)/2);
    C = 2 * gamma(d/2) / (sqrt(pi) * gamma((d-1)/2));
    val = C * 1/x * integral(ff,0,x);
  
  function val = myp(t)
    val = exp(-(t.^2)./5);
  
  function val = hbar(x,t)
    val = exp(myp(t) .* cos(2*pi*x) .* cos(myp(t).*sin(2*pi*x)))./exp(myp(0));
    

