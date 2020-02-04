function [x,crit] = adaptdesc1D(x,f,g, fft_f, t, freq, N, n)

maxiter = 5000;
G = feval(g,x, fft_f, t, freq, N, n);
stepsize = 1/norm(G(:));
stepmult = 1.2;
stepdiv = 2;
breakratio = 1e-10;
verbosemode = 0;

crit(1) = feval(f,x, fft_f, t, freq, N, n);
for iter = 1:maxiter
    G = feval(g,x,  fft_f, t, freq, N, n);
    stepsize = stepsize * stepdiv;
    minimtest = 0;
    while(~minimtest) 
        stepsize = stepsize / stepdiv;
        xnew = x - stepsize * G'; 
        xnew(1) = -sum(xnew(2:n)); % Condition d'identifiabilite
        crit(iter+1) = feval(f,xnew, fft_f, t, freq, N, n);
        minimtest = (crit(iter+1) < crit(iter));
    end
    x = xnew;
    
    stepsize = stepsize * stepmult;
    if verbosemode
        disp(['iteration ',num2str(iter),...
            '     functional = ',num2str(crit(iter)),...
            '     stepsize = ', num2str(stepsize)]);
    end
    % termination criteria
    if (crit(iter)-crit(iter+1)) < breakratio*(crit(1)-crit(iter+1))
        sum(x)
        return
    end
end

if verbosemode
    disp('maximum number of iterations exceeded');
end
