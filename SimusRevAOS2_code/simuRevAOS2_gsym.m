% Simulations shift aleatoires
addpath('/Users/bgemily/Documents/Academic/SC/WavDen/codes')
addpath('/Volumes/jbigot/Travail_iDisk/Programmation/SimusPaperRandomShifts/SimusRevAOS2_code')
cd('/Volumes/jbigot/Travail_iDisk/Programmation/SimusPaperRandomShifts/SimusRevAOS2_code')

clear all;
close all;

N = 256;
t = linspace(0,1,N);


namefun = 'Bumps';

if strcmp(namefun,'Gauss')
    fstar = 0.7*exp(-(t-0.4).^2/(2*0.025^2)) + 0.4*exp(-(t-0.2).^2/(2*0.01^2)) + exp(-(t-0.7).^2/(2*0.045^2));
elseif strcmp(namefun,'Sine')
    fstar = MakeSignalNewb('Time Shifted Sine',N);
elseif strcmp(namefun,'Wave')
    fstar = MakeSignalNewb('Wave',N);
elseif strcmp(namefun,'Bumps')
    fstar = MakeSignalNewb('Bumps',N);
elseif strcmp(namefun,'HeaviSine')
    fstar = MakeSignalNewb('HeaviSine',N);
elseif strcmp(namefun,'Blocks')
    fstar = MakeSignalNewb('Blocks',N);
elseif strcmp(namefun,'Doppler')
    fstar = MakeSignalNewb('Doppler',N);
elseif strcmp(namefun,'Creneau')
    fstar=base(t);
else
    error('Mauvais nom pour namefun')
end


figure(1)
plot(t,fstar,'LineWidth',2)
% eval(['print -depsc ' namefun '.eps ;'])

n = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Choix du type de bruit
%%%%%%%%%%%%%%%%%%%%%%%%%

% % Gaussien
% sig = 0.035;
% theta = sig*randn(1,n);
% truedens = normpdf(t,0,sig)+normpdf(reverse(t),0,sig);
% 
% % Gaussian mixture
% mu1 = -0.3;
% mu2 = 0.1;
% sig1 = 0.03;
% sig2 = 0.03;
% p1 = 1/4;
% p2 = 1-p1;
% sigmax = sqrt(p1*(sig1^2+mu1^2) + p2*(sig2^2+mu2^2) - (p1*mu1 + p2*mu2)^2);
% c = binornd(1,p1,1,n);
% theta = c.*normrnd(mu1,sig1,1,n) + (1-c).*normrnd(mu2,sig2,1,n);
% truedens = p1*normpdf(t,mu1,sig1) + p2*normpdf(t,mu2,sig2) + p1*normpdf(reverse(t),mu1,sig1) + p2*normpdf(reverse(t),mu2,sig2) ;

% % Vraies valeurs de gamma
% omega = 2*pi*([0:(N/2) (-N/2+1):-1]);
% gamma = p1*exp(-mu1*i*omega-sig1^2*omega.^2/2) + p2*exp(-mu2*i*omega-sig2^2*omega.^2/2);
% 
% figure(8)
% plot(1:N, abs(fft(fstar)), 1:N, abs(fft(fconv)./gamma),'--r')

% Laplace
mu = 0.03;
theta = exprnd(mu,1,n)-exprnd(mu,1,n);
m = 0;
truedens = (exppdf(t-m,mu)+exppdf(flip(t-m),mu))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the noisy curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnoise = zeros(n,N);
rsnr = 7;
s = std(fstar)/rsnr;

% On utilise des variances différentes pour chaque courbe
b = 0.09*s;
sigma = s+-b+2*b*rand(1,n);

for (m = 1:n)
    fnoise(m,:) = decale(fstar,theta(m)) + sigma(m)*randn(1,N);
end

%  Visualisation des courbes decalees
figure(2)
plot(t,fnoise(1:10,:),'LineWidth',2)
namefig = sprintf('samplen%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

% Moyenne naive
naivemean = mean(fnoise);
figure(3)
plot(t,fstar,'--b',t,naivemean,'r','LineWidth',2)
namefig = sprintf('naiven%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

% Calcul de la fft des courbes bruitées
fft_fnoise = fft(fnoise,[],2);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation de sigma
%%%%%%%%%%%%%%%%%%%%%%%%%

hat_sigma = zeros(1,n);
qmf = MakeONFilter('Symmlet',8);

for (m = 1:n)
    [y,coef] = NormNoise(fnoise(m,:),qmf);
    hat_sigma(m) = 1/coef;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation des shifts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Fourier transform for frequency lower than freq_max = l0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Smoothing parameter
freq_max = 3;
freq = (2*pi) * ((-freq_max):freq_max); % frequency along the x-axis

fft_fnoise_freq_max = zeros(N,length(freq));

for (m = 1:n)
    fft_fnoise_freq_max(m,:) = fourier1D(fnoise(m,:), t, freq); %%% Manual Fourier transformation
end

%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation
%%%%%%%%%%%%%%%%%%%%%%%

% Gradient descent in the fourier domain with denoising
hat_theta_init = zeros(1,n);

%%% F1D: objective function, the mean squared distance between shifted curves and
% their center
[hat_theta ,crit] = adaptdesc1D(hat_theta_init,'F1D','GradF1D', fft_fnoise_freq_max, t, freq, N, n);
hat_theta = -hat_theta;

u = linspace(min(theta),max(theta),1000);

figure(4)
plot(hat_theta,theta,'r*',u,u,'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation de gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 2*pi*([0:(N/2) (-N/2+1):-1]);
gamma = mean(exp(-i*theta(2:n)'*omega));
hat_gamma = mean(exp(-i*hat_theta(2:n)'*omega));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation par ondelettes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = floor(log2(N));
j0 = 3;

% Parameters for the Meyer wavelets
deg = 3;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the wavelet coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deconvolution step
Nfreq = 40;
weight = 1./hat_gamma;
weight((Nfreq+1):(N-Nfreq)) = 0;
hatf_fft1 = mean(fft_fnoise).*weight;
f_fft = fft(fstar);

% Estimation par recalage
hatf_fft2 = mean(fft_fnoise.*exp(i*hat_theta'*omega));

figure(5)
plot(t,fstar,t,real(ifft(hatf_fft1)),'r',t,real(ifft(hatf_fft2)),'c')

%  Compute Coefficients at Coarse Level
wtrue = zeros(1,N);
wnoise1 = zeros(1,N);
wnoise2 = zeros(1,N);
wtrue(1:(2^j0)) = CoarseMeyerCoeff(f_fft,j0,N,deg);

% Estimation par deconvolution
wnoise1(1:(2^j0)) = CoarseMeyerCoeff(hatf_fft1,j0,N,deg);

% Estimation par recalage
wnoise2(1:(2^j0)) = CoarseMeyerCoeff(hatf_fft2,j0,N,deg);

%  Loop to Get Detail Coefficients for levels  j=j0,...,J-2.
for j = j0:(J-2),
  wtrue1((2^j+1):(2^(j+1))) = DetailMeyerCoeff(f_fft,j,N,deg);
  wnoise1((2^j+1):(2^(j+1))) = DetailMeyerCoeff(hatf_fft1,j,N,deg);
  wnoise2((2^j+1):(2^(j+1))) = DetailMeyerCoeff(hatf_fft2,j,N,deg);
end

%  Calculate Fine Level Detail Coefficients (for j=J-1).
wtrue((2^(J-1)+1):(2^J)) = FineMeyerCoeff(f_fft,N,deg);
wnoise1((2^(J-1)+1):(2^J)) = FineMeyerCoeff(hatf_fft1,N,deg);
wnoise2((2^(J-1)+1):(2^J)) = FineMeyerCoeff(hatf_fft2,N,deg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sigma_jk
Sigma_jk = zeros(1,N);

for (j = j0:(J-1))
    pos = dyad(j);
    for (k = 1:2^j)
        wMeyer = zeros(1,N);
        wMeyer(pos(k)) = 1;
        psi = IWT_YM(wMeyer,j0,deg);
        fft_psi = fft(psi)/N^2;
        Sigma_jk(pos(k)) = sum( (abs(fft_psi)./abs(hat_gamma)).^2 );
    end
end

% Threshold
eta = 1;
sig = sqrt(mean(hat_sigma.^2));
thr1 = 2*sig*sqrt(Sigma_jk) * sqrt(2*eta*log(n)/n);
thr2 = 2*sig*sqrt(Sigma_jk) * sqrt(2*eta*log(n)/n);

figure(6)
indice = (2^j0+1):N;
plot(indice,abs(wtrue(indice)),'--b',indice,abs(wnoise1(indice)),'r',indice,thr1(indice),'r',indice,abs(wnoise2(indice)),'c',indice,thr2(indice),'c','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thresholding of the wavelet coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hatw1 = wnoise1;
hatw2 = wnoise2;
thr_type = 'H'; % Type of thresholding either Hard or Soft
j1 = J-1;

if strcmp(thr_type,'S')
    for j=j0:j1,
        pos = dyad(j);
        for (k = 1:2^j)
            hatw1(pos(k)) = SoftThresh(wnoise1(pos(k)),thr1(pos(k))) ;
            hatw2(pos(k)) = SoftThresh(wnoise2(pos(k)),thr2(pos(k))) ;
        end
    end
    hatw((2^(j1+1)+1):(2^J)) = 0;
else
     for j=j0:j1,
        pos = dyad(j);
        for (k = 1:2^j)
            hatw1(pos(k)) = HardThresh(wnoise1(pos(k)),thr1(pos(k))) ;
            hatw2(pos(k)) = HardThresh(wnoise2(pos(k)),thr2(pos(k))) ;
        end
    end
    hatw1((2^(j1+1)+1):(2^J)) = 0;
    hatw2((2^(j1+1)+1):(2^J)) = 0;
end

% Estimation with wavelet thresholding
hatf1 = IWT_YM(hatw1,j0,deg);
hatf2 = IWT_YM(hatw2,j0,deg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation avec algo min-min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moyenne naive
naivemean = mean(fnoise);

% Nombre d'iterations de notre algorithme
iter=3;
Z = naivemean;
hat_theta = zeros(size(theta));

% Boucle d'iterations de l'algorithme
for i=1:iter
    i
    % Boucle sur le nombre de courbes
    aux=zeros(1,length(Z));
    for m=1:n
       hat_theta(m)= eval_shift(Z,fnoise(m,:));
       aux = aux + decale(fnoise(m,:),-hat_theta(m));
    end
    Z=aux/n;
end
hatf_minmin = Z;

figure(8)
plot(t,fstar,'--b',t,hatf1,'r','LineWidth',2)
namefig = sprintf('hatf1%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

figure(9)
plot(t,fstar,'--b',t,hatf2,'r','LineWidth',2)
namefig = sprintf('hatf2%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

figure(10)
plot(t,fstar,'--b',t,hatf_minmin,'r','LineWidth',2)
namefig = sprintf('moyminminn%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

