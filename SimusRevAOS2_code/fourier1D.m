function y = fourier1D(f, t, freq)

clear('i')
N = length(t);
y = (exp(-i.*(freq'*t)) * f')' ./ N;

