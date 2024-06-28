clc
clearvars
close all
%--------------------------------------------------------------------------
%---------------------------input------------------------------------------
G=input('(binery=1, binery phase=2, sinusidal=3, sinusidal phase=4):');%Gtating
B=input('(Gaussian Beam=1,Besel Beam=2,Lagger-Gaussian Beam=3,airy Beam=4,airy Beam=5):');%Beams
P=input('(ASM=1,HCM=2):');%Angular Spectrum or Hugence Convolution
lxx =input('lx,periodic for line in grating(micronmetre(1e-6)) :');
zz=input('z,Leangh propagation (m):');%Leangh propagation  of z
%------------------------------Real Space----------------------------------
Nx = 2^10;
Ny = 2^10;
lx = lxx*(1e-6);

x0_min = -3e-3;
x0_max =  3e-3;
y0_min = -3e-3;
y0_max =  3e-3;
x0 = linspace(x0_min,x0_max,Nx);
y0 = linspace(y0_min,y0_max,Ny);

dxo = x0(2) - x0(1);
dyo = y0(2) - y0(1);

[x,y] = meshgrid(x0,y0);
r = sqrt(x.^2+y.^2);
rho = atan(y./x);

L_D=2*(abs(x0_max)+abs(x0_min))/lx;%Number of dark and light
PP=Nx/(L_D);%pixl of dark and light bars
length=PP*dxo;%length of dark and light bars
total=length*PP;%total length of gratting,
%The number and size of dark and light are one because our gratting is a sin

%-----------------------------Input Field----------------------------------
lambda = 632.8e-9;
k = 2.*pi/lambda;

% ---------------------------Gratings--------------------------------------
% -------- Sin
if G==1
    % -------- binery
    Grating = sin((2*pi*x/lx));
    Grating = heaviside(Grating);
    %----------
elseif G==2
    %-------binery phase
    Grating = sin((2*pi*x/lx));
    Grating = heaviside(Grating);
    Grating = exp(1i.*Grating);
    
    %----------
elseif G==3
    %------------- sinusidal  grating
    Grating = 0.5*(sin((2*pi*x/lx))+1);
    %----------
elseif G==4
    %------------------- sinusidal phase grating
    Grating = 0.5*(sin((2*pi*x/lx))+1);
    Grating = exp(1i.*Grating);
    %----------
end
% ylim([-5 5]);
% xlim([0 15]);

% ---------------------------Apartures-------------------------------------
Ap = zeros(size(r));
% -------- Circular
Ap(r <= 2000e-6) = 1;
% ----------------------------Beams----------------------------------------
if B==1
    % -------- Gaussian Beam
    wg = 700e-6;
    UU= exp(-(r/wg).^2);%U_GB
    U= mnormalize(UU);
    %----------
    
elseif B==2
    %----------------Besel Beam
    w_B1 = 50e-6;w_B2 = 1e-3;
    UU= besselj(0,(r/w_B1)).*exp(-(r/w_B2).^2);%U_BB
    U= mnormalize(UU);
    %----------
elseif B==3
    %-----------------Lagger-Gausian Beam
    p = 1;                  % Degree
    l = 2;                  % Order
    w_LG0 = 50e-6;             % Beam waist
    z_R = k*w_LG0^2/2;      % Calculate the Rayleigh range
    z_LG = 0.0;
    [phi_LG, r_LG] = cart2pol(x, y);    %Transform Cartesian coordinates to polar
    U00 = 1/(1 + 1i*z_LG/z_R) .* exp(-r_LG.^2/w_LG0^2./(1 + 1i*z_LG/z_R));
    w_LG = w_LG0 * sqrt(1 + z_LG.^2/z_R^2);
    R_LG = sqrt(2)*r_LG./w_LG;
    Lpl = nchoosek(p+l,p) * ones(size(R_LG));   % x = R(r, z).^2
    % com = Binomial coefficient or all combinations
    for m = 1:p
        Lpl = Lpl + (-1)^m/factorial(m) * nchoosek(p+l,p-m) * R_LG.^(2*m);
    end
    
    UU = U00.*R_LG.^l.*Lpl.*exp(1i*l*phi_LG) ...
        .*exp(-1i*(2*p + l + 1)*atan(z_LG/z_R));%U_LG
    U= mnormalize(UU);
    %----------
elseif B==4
    
    %-------------------------------------1D Airy Beam
    a = 0.25;    % truncation coefficient
    s = .18;     % mainlop width (mm)
    
    D = (airy(x./s)).*exp(a.*x);   % truncated airy beam
    
    BB = D.*conj(D);
    
    up=BB-min(BB(:));
    UU=up./max(up(:));
    U= mnormalize(UU);
    %----------
elseif B==5
    
    %------------------------------2D Airy Beam
    a = 0.05;    % truncation coefficient
    s = .15;     % radius
    r0 = 0.4;    % mainlop width (mm)
    
    rr = (x.^2+y.^2);
    D  = airy((r0 - rr)./s).*exp(a.*((r0-rr)./s));
    
    BB = D.*conj(D);
    up=BB-min(BB(:));
    UU=up./max(up(:));
    U= mnormalize(UU);
    %----------
end
%--------------------------------Wavefront---------------------------------
% --------------------------Plane
phi = zeros(size(r));

%---------------------------------Beam Interaction  by Grating-------------

U00 = U.*Ap.*  Grating ;
U0= mnormalize(U00);
I0 = abs(U0).^2;
figure(1);
imagesc(x0*1e3, y0*1e3, I0)

xlabel('$x\ \textrm{[cm]}$','interpreter','latex','FontSize',22);
ylabel('$y\ \textrm{[mm]}$','interpreter','latex','FontSize',22);
ylabel(colorbar,'$ \textrm{Intensity\ [arb. u.]}$','FontSize',22,...
    'interpreter','latex');
axis image;
colorbar;

%-----------------------------Propagation----------------------------------

Ini.dxo = dxo;
Ini.dyo = dyo;
Ini.lambda = lambda;
Nz = 50;
z = linspace(0,zz,Nz);

for n = 1 : Nz
    d = z(n);
    
    if P==1
        
        Ud_ASM = Propagate_ASM(U0,d,Ini);
        %The Angular Spectrum method is used for short  distances
        %(such as Fresnel  diffraction which is used for short  distances.)
        
        I_ASM = abs(Ud_ASM).^2;
        II(:,n) = I_ASM(round(end/2),:);%Intensity
        imagesc(x0*1e3, y0*1e3, I_ASM);
    elseif P==2
        
        Ud_HCM = Propagate_HCM(U0,d,Ini);
        %The Hugence Convolution  method is used for long distances
        %(such as Fraunhofer  diffraction which is used for long distances.)
        
        I_HCM = abs(Ud_HCM).^2;
        II(:,n) = I_HCM(round(end/2),:);%Intensity
        imagesc(x0*1e3, y0*1e3, I_HCM);
    end
    
    
    title(fix(100*d/max(z))*max(z)/100)
    xlabel('$x\ \textrm{[mm]}$','interpreter','latex','FontSize',22);
    ylabel('$y\ \textrm{[mm]}$','interpreter','latex','FontSize',22);
    axis image;
    colorbar;
    
    pause(.001)
end

figure;
hold on
%-----------Plotting amplitude at the initial and final planes-------------
if P==1
    plot(y0*1e3,I_ASM(:,round(end/2)))
elseif P==2
    plot(y0*1e3,I_HCM(:,round(end/2)))
end

xlabel('$x \textrm{[mm]}$','interpreter','latex','FontSize',22);
ylabel('$I \textrm{[arb.\ u.]}$','interpreter','latex','FontSize',22);
grid on;

%----------------------------------y-z Intensity---------------------------
figure;

imagesc(z*1e3,y0*1e3,II)

xlabel('$z\ \textrm{[mm]}$','interpreter','latex','FontSize',22);
ylabel('$y\ \textrm{[mm]}$','interpreter','latex','FontSize',22);
ylabel(colorbar,'$ \textrm{Intensity\ [arb. u.]}$','FontSize',22,...
    'interpreter','latex');

%--------------------------------------------------------------------------