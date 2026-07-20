function [alpha, x] = resolve(A,b)
    [m,n]=size(A);
    size(A);
    Im = sparse(speye(m));
    In = sparse(speye(n));
    M = sparse(m+n,m+n);
    M = [Im,A;-A',In];
    nb = zeros(m+n,1);
    nb(1:m) = b;
    x = M\nb;
    alpha = x(m+1:end);
end
