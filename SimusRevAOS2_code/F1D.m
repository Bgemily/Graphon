function y = F1D(theta, fft_f, t, freq, N, n)

clear('i')

fft_f_tilde = zeros(n,length(freq));

for (m = 1:n)
    fft_f_tilde(m,:) = fft_f(m,:) .* exp(i*theta(m)*freq);
end

totmoy = 0;
for (m = 1:n)
    totmoy = totmoy + fft_f_tilde(m,:);
end
totmoy = totmoy/n;

tot = zeros(1,n);
for (m = 1:n)
    tot(m) = sum(sum( abs(fft_f_tilde(m,:) - totmoy).^2 ));
end
y = sum(tot)/n;