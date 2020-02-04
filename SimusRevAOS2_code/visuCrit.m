% Simulations shift aleatoires
addpath('/users/jeremiebigot/Travail/Programmation/WavDen/codes')
addpath('/Volumes/jbigot/Travail_iDisk/Programmation/SimusPaperRandomShifts/SimusRevAOS2_code')
cd('/Volumes/jbigot/Travail_iDisk/Programmation/SimusPaperRandomShifts/SimusRevAOS2_code')

clear all;
close all;

N = 256;
t = linspace(0,1,N);


namefun = 'Creneau';

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

n = 2;

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
mu = 0.05;
theta = exprnd(mu,1,n)-exprnd(mu,1,n);
m = 0;
truedens = (exppdf(t-m,mu)+exppdf(reverse(t-m),mu))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the noisy curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnoise = zeros(n,N);
rsnr = 7;
s = std(fstar)/rsnr;

% On utilise des variances différentes pour chaque courbe
b = 0.09*s;
sigma = s+-b+2*b*rand(1,n);

sigma = zeros(1,n);

for (m = 1:n)
    fnoise(m,:) = decale(fstar,theta(m)) + sigma(m)*randn(1,N);
end

%  Visualisation des courbes decalees
figure(2)
plot(t,fnoise(1:2,:),'LineWidth',2)
namefig = sprintf('samplen%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

% Moyenne naive
naivemean = mean(fnoise);
figure(3)
plot(t,fstar,'--b',t,naivemean,'r','LineWidth',2)
namefig = sprintf('naiven%d.eps',n);
% eval(['print -depsc ' namefun namefig ';'])

% Smoothing parameter
freq_max = floor(N/2);
freq = (2*pi) * ((-freq_max):freq_max); % frequency along the x-axis

fft_fnoise_freq_max = zeros(N,length(freq));

for (m = 1:n)
    fft_fnoise_freq_max(m,:) = fourier1D(fnoise(m,:), t, freq);
end


M = 100;
u = linspace(-0.25,0.25,M);
D = zeros(1,M);

for (k=1:M)
    D(k) = F1D([u(k) -u(k)], fft_fnoise_freq_max, t, freq, N, n);
end

figure(1)
thetastar = theta-mean(theta);
plot(u,D,thetastar,zeros(1,n),'r*','LineWidth',2)
