function X = rand_power_law_gaussian(d,ni,nj,ymax,ymin,ymean,ystd)

% Random numbers y (matrix ni x nj) between ymin and ymax
% with a probability density p(y) proportional to |y|^d * p_normal(y,ymean,ystd)
%
% Example: hist(rand_power_law_gaussian(0,1000,1,2,0,1,0.1),100);

eval(default('d','1','ni','1','nj','1','ymax','1','ymin','0','ymean','0','ystd','1'));
if ymax < ymin, error('ymax < ymin; reverse order of arguments'); end

pdf_function = inline(sprintf('exp(-(x-%f)^2/%f^2) * abs(x)^%d',ymean,ystd,d));
x = rand_given_pdf(pdf_function, ymin, ymax, ni, nj)