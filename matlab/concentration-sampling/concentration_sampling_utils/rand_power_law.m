function X = rand_power_law(d,ni,nj,ymax,ymin)

% Random numbers y (matrix ni x nj) between ymin and ymax
% with a probability density p(y) proportional to |y|^d
%
% Example: hist(rand_power_law(3,10000,1,2,-1),100);

eval(default('d','1','ni','1','nj','1','ymax','1','ymin','0'));

if [ymin == 0] * [ymax >= 0],
  X = ymax * [rand(ni,nj).^[1/[d+1]]];
elseif [ymax == 0] * [ymin <= 0],
  X = ymin * [rand(ni,nj).^[1/[d+1]]];
else
  ym = max(abs(ymin),abs(ymax));
  X = ym * sign(rand(ni,nj)-0.5) .* [rand(ni,nj).^[1/[d+1]]];
  ind_replace = find( [X < ymin] + [X > ymax] );
  nr = length(ind_replace);
  while nr,
    X(ind_replace) = ym * sign(rand(nr,1)-0.5) .* [rand(nr,1).^[1/[d+1]]];
    ind_replace = find( [X < ymin] + [X > ymax] );
    nr = length(ind_replace);
  end
end
