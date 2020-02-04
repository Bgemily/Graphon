function res = GradF1D(theta, fft_f, t, freq, N, n)

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

res = zeros(n,1);
for (m = 1:n)
    res(m) = -(2/n)* sum( real( (i * freq).* fft_f_tilde(m,:).*conj(totmoy) ));
end

