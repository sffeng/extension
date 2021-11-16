function plot_turning
  rng(1985)
  close all
  n = 60; %n=20 and rnd1985 pretty good
  d = 2;

  figure

  subplot(2,3,1)
  t=0;
  plothdxt(n,t,d)
  hold on
  plotunitcircle()
  axis square
  ylabel('Covariance')

  subplot(2,3,4)
  [X,Y,S] = getcorrS(n,t,d);
  mesh = [X(:) Y(:)];
  corr.name = 'gauss';
  corr.c0 = [1];
  corr.sigma=S(:);
  [F,KL] = randomfield(corr,mesh,'trunc',10);
  surf(X,Y,reshape(F,n,n))
  hold on
  plotunitcircle()  
  view(2)
  shading interp
  hold on
  plotunitcircle()
  axis square
  ylabel('Realization')
  

  subplot(2,3,2)
  t=1;
  plothdxt(n,t,d)
  hold on
  plotunitcircle()
  axis square

  subplot(2,3,5)
  [X,Y,S] = getcorrS(n,t,d);
  mesh = [X(:) Y(:)];
  corr.name = 'gauss';
  corr.c0 = [1];
  corr.sigma=S(:);
  [F,KL] = randomfield(corr,mesh,'trunc',10);
  surf(X,Y,reshape(F,n,n))
  hold on
  plotunitcircle()  
  view(2)
  shading interp
  hold on
  plotunitcircle()
  axis square
%  colorbar


  subplot(2,3,3)
  t=2;
  plothdxt(n,t,d)
  hold on
  plotunitcircle()
  axis square
  
  subplot(2,3,6)
  [X,Y,S] = getcorrS(n,t,d);
  mesh = [X(:) Y(:)];
  corr.name = 'gauss';
  corr.c0 = [1];
  corr.sigma=S(:);
  [F,KL] = randomfield(corr,mesh,'trunc',10);
  surf(X,Y,reshape(F,n,n))
  hold on
  plotunitcircle()  
  view(2)
  shading interp
  hold on
  plotunitcircle()
  axis square
  
  
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
    

