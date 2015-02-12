function f = ecm_regularisation(x,x_min,x_max,lambda)

f = lambda * 0.5 * sum([x - 0.5 * [x_min + x_max]].^2);
