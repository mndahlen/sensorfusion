a = 1;
b = -2;
c = 3;
N = 1000;
y = a + b*log((1:N)' + c);

cg = 0:0.1:5;
ahat = zeros(size(cg));
bhat = zeros(size(cg));
V = zeros(size(cg));
for i = 1:length(cg)
    H = [ones(N,1) log((1:N)' + cg(i))];
    xhat = H\y;
    ahat(i) = xhat(1);
    bhat(i) = xhat(2);
    yhat = ahat(i) + bhat(i).*log((1:N)' + cg(i));
    V(i) = 0.5*sum((y-yhat).^2);
end
[~, ind] = min(V);
ahat(ind)
bhat(ind)
cg(ind)

