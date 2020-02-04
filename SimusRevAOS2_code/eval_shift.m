function u=eval_shift(f,g)

% Via fminsearch
u = fminsearch(@(theta) crit(theta,f,g),0);



