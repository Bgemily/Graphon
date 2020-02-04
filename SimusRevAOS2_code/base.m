function res=base(x)
%Fonction creneau de base
res = zeros(size(x));
res((x>=2/5) & (x<=2/5+0.2)) = 1;