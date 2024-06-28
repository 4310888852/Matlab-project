function U = Propagate_ASM(U0,Dz,Ini) 
%The Angular Spectrum Method to propagate waves to diffraction near field.  

% Initial values of incident wave 
dxo = Ini.dxo; % The pixel size of plane along x axis
dyo = Ini.dyo; % The pixel size of plane along y axis
lambda = Ini.lambda; % The wavelength of incident wave

[Ny,Nx] = size(U0); %The number of pixels 
Lx = Nx*dxo; % The length of plane along x axis
Ly = Ny*dyo; % The length of plane along y axis

% The coordinates of Fourier domain
fx = (-1/(2*dxo) : 1/Lx : 1/(2*dxo)-1/Lx)-mod(Nx,2)/Lx/2; 
fy = (-1/(2*dyo) : 1/Ly : 1/(2*dyo)-1/Ly)-mod(Ny,2)/Ly/2; 
[fx,fy] = meshgrid(fx,fy);

% H = exp(-1j*pi*lambda*Dz*(fx.^2+fy.^2)); %transfer function
H = exp(1j*2*pi*Dz*sqrt(1/lambda.^2-fx.^2-fy.^2)); %transfer function
H = fftshift(H); 
U1 = fft2(U0); 
U2 = H.*U1; 
U = ifft2(U2);
