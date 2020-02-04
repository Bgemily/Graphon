function cost = crit(theta,f,g)

    fdecal = decale(f,theta);   
    cost =sum((fdecal-g).^2);