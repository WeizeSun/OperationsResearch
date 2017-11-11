function [f,df,ddf] = obj(x)

n = length(x);

y    = (x-1)/4;
y_sq = y.^2;
s    = sin(pi*(1+y));
s_sq = s.^2;
c    = cos(pi*(1+y));

f = pi/n * (10*s_sq(1) + sum(y_sq(1:n-1) .* (1 + 10*s_sq(2:n))) + y_sq(n));

if nargout > 1
  dg_dy = 2*y .* [1+10*s_sq(2:n); 1];
  dg_ds = 20*s .* [1; y_sq(1:n-1)];
  ds_dx = pi/4 * c;

  df = pi/n * (dg_dy/4 + dg_ds .* ds_dx);
  
  if nargout > 2
    ddg_dydy = diag([2+20*s_sq(2:n); 2]);
    ddg_dsds = diag([20; 20*y_sq(1:n-1)]);
    ddg_dsdy = diag(40*y(1:n-1).*s(2:n), -1);
    ddg_dyds = ddg_dsdy.';
    dds_dxdx = diag(-pi^2/16*s);
    
    ddf = pi/n * ...
        ((ddg_dydy/4 + ddg_dsdy .* repmat(ds_dx,1,n))/4 + ...
         (ddg_dyds/4 + ddg_dsds .* repmat(ds_dx,1,n)) .* repmat(ds_dx,1,n).' + ...
         repmat(dg_ds,1,n).' .* dds_dxdx);
  end
end