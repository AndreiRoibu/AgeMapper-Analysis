function [coef, p, n, z] = nancorr(A, B)
  %NANCORR - Pearson correlation coefficient
  %
  % B can be left out
  %
  %  coef = NANCORR(A, B) is equivalent to 
  %  coef = corr(A, B, 'rows','pairwise'),
  %  but NANCORR works much faster.
  %  
  %  [coef, p, n, z] = NANCORR(A, B)
  %  INPUT:
  %    A, B - input matrices, single or double, with equal number of rows
  %
  %  OUTPUT:
  %    coef - matrix of Pearson correlation coefficients
  %    p    - matrix of p-values
  %    n    - matrix containing the number of defined values
  %    z    - matrix of fisher-transformed correlations
  %
  % NOTES  
  %  pvalue can be calculated as 2*tcdf(-abs(t), n - 2)

  %Am, Ap: m=missing elements, p=present elements

  if ( exist('B') == 0 )
    B=A;
  end

  Am=~isfinite(A); Bm=~isfinite(B); 
  
  if strcmp(class(A), 'single')
    Ap=single(~Am); Bp=single(~Bm);    
  else
    Ap=double(~Am); Bp=double(~Bm);
  end

  % zero out nan elements
  A(Am)=0; B(Bm)=0;

  % code one of the formulas from https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
  % this procedure might be numericaly unstable for large values,
  % it might be reasonable to center each column before calling nancorr.

  xy = A' * B;          % sum x_i y_i
  n  = Ap' * Bp;        % number of items defined both in x and y
  mx = A' * Bp ./ n;    % mean values in x, calculated across items defined both in x and y
  my = Ap' * B ./ n;    % mean values in y, calculated across items defined both in x and y
  
  x2 = (A.*A)' * Bp;    % sum x^2_i, calculated across items defined both in x and y
  y2 = Ap' * (B.*B);    % sum y^2_i, calculated across items defined both in x and y
  
  % sx, sy - standard deviations 
  % sx   = sqrt(x2 - n .* (mx.^2));
  % sy   = sqrt(y2 - n .* (my.^2));
  sx = x2 - n .* (mx.^2); sx(sx<0)=0; sx=sqrt(sx);
  sy = y2 - n .* (my.^2); sy(sy<0)=0; sy=sqrt(sy);
  
  coef = (xy - n .* mx .* my) ./ (sx .* sy);      % correlation coefficient
  coef(n<3)=nan;

  t    = coef .* sqrt((n - 2) ./ (1 - coef.^2));  % t-test statistic
  p    = 2*tcdf(-abs(t), n - 2);   p(p==0)=1e-300;

  grot1=1+coef; grot1(grot1<=0)=1e-100; % handle cases where r<=-1 or >=1
  grot2=1-coef; grot2(grot2<=0)=1e-100;
  z = sqrt(max(n-3,1)) .* (0.5*log(grot1./grot2));

end

