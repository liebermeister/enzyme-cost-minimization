% ----------------------------------------------------------------------------------
% plot enzyme costs of a simple mass action kinetics
% ----------------------------------------------------------------------------------

clear
ca

v_list = -2:0.1:-0.1;
nc = length(v_list);

for it1 = 1:nc,
  for it2 = 1:nc
    x(it1,it2) = v_list(it1);
    y(it1,it2) = v_list(it2);
  end
end

a = 1+exp(1/2 * [x+y]) ./ [1-exp(1/2*[x+y])];
b = 1/2 * [ [1+exp(x)]./[1-exp(x)] + [1+exp(y)]./[1-exp(y)]];
im(a<b)