function U = Propagate_HCM(U0,Dz,Ini)
%The Convolution Huygens Method to propagate waves to diffraction far field.  

% Initial values of incident wave 
dxo = Ini.dxo; % The pixel size of plane along x axis
dyo = Ini.dyo; % The pixel size of plane along y axis
lambda = Ini.lambda; % The wavelength of incident wave
k = 2*pi/lambda; 

[Ny,Nx] = size(U0); %The number of pixels 
Lx = Nx*dxo; % The length of plane
Ly = Ny*dyo; % The length of plane

x = -Lx/2+dxo : dxo : Lx/2;
y = -Ly/2+dyo : dyo : Ly/2;
[x,y] = meshgrid(x,y);

h = 1/(1i*lambda*Dz) .* exp(1i*k*sqrt(x.^2+y.^2+Dz.^2)); %Impuls Respons
H = fft2(fftshift(h)) * dxo^2; 
U1 = fft2(fftshift(U0)); 
U2 = H .* U1; 
U = ifftshift(ifft2(U2));

