function res = decale(f,theta)

    N = length(f);
    fext = [f f f];
    u = linspace(-1,2,3*N);
    fdecal = interp1(u,fext,u-theta,'linear',0);
    res = fdecal((N+1):(2*N));