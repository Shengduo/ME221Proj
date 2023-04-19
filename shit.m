%  Use of 2D FFT for Fourier-spectral solution of
%
%      u_xx + u_yy  = (r^2 - 2*sig^2)*exp(-r^2/(2*sig^2))/sig^4,  0 < x,y < 1
%
%  where
%
%      r^2 = (x-0.5)^2 + (y-0.5)^2, sig = 0.2
%
%  with periodic BCs on u in x,y, using N = 16 modes in each direction.
%  For sig << 1, this Poisson eqn has exact solution
%
%      uex = exp(-r^2/(2*sig^2)) + C (C is an arbitrary constant)
%
%  We will choose C by taking u(x=0,y=0) = 0. In the limit sig << 1, this
%  implies C = 0.
%  Script makes a surface plot of u at the Fourier grid points.
  clc,clear;
  close all;
  N = 2048;      % No. of Fourier modes...should be a power of 2
  L = 1;       % Domain size (assumed square)
  sig = 0.1;   % Characteristic width of f (make << 1)

  k = (2*pi/L)*[0:(N/2-1) (-N/2):(-1)]; % Vector of wavenumbers
  [KX KY]  = meshgrid(k,k); % Matrix of (x,y) wavenumbers corresponding
                            % to Fourier mode (m,n)
  delsq = -(KX.^2 + KY.^2); % Laplacian matrix acting on the wavenumbers
  delsq(1,1) = 1;           % Kluge to avoid division by zero of (0,0) waveno.
                            % of fhat (this waveno. should be zero anyway!)
  
  % Construct RHS f(x,y) at the Fourier gridpoints

  h = L/N;     % Grid spacing
  x = (0:(N-1))*h ;
  y = (0:(N-1))*h;
  [X, Y] = meshgrid(x,y);
  rsq = (X-0.5*L).^2 + (Y-0.5*L).^2;
  sigsq = sig^2;
  f = exp(-rsq/(2*sigsq)).*(rsq - 2*sigsq)/(sigsq^2);

%  Spectral inversion of Laplacian

  fhat = fft2(f);
  u = real(ifft2(fhat./delsq));
  u = u - u(1,1);   % Specify arbitrary constant by forcing corner u = 0.
  
  lapalacian_u = (u([2:N, 1], :) - 2 .* u + u([N, 1:N-1], :)) / h^2 ...
                 + (u(:, [2:N, 1]) - 2 .* u + u(:, [N, 1:N-1])) / h^2;
  residual = sum((lapalacian_u - f) .^ 2 * h^2, 'all');
  
%  Plot out solution in interior
  disp(strcat("residual = ", num2str(residual)));
  uex = exp(-rsq/(2*sigsq));
  errmax = norm(u(:)-uex(:),inf)
%   surf(X,Y,u)
%   xlabel('x')
%   ylabel('y')
%   zlabel('u')
%   title('Fourier spectral method for 2D Poisson Eqn')