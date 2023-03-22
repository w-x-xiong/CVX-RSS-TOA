function [x_est, fail] = HuberUnder(d0, r, P0, P, gamma, anc, R)

fail = false;

[H, L] = size(anc);

x_est = zeros(H,1);
beta = zeros(L,1);

for i = 1:L
    beta(i) = 10^((P(i)-P0)/10/gamma);
end

cvx_begin quiet
variables x_var(H)
expression obj
for i = 1:L
    obj = obj + huber_pos((10*gamma*beta(i)*norm(x_var - anc(:,i))/d0/log(10) - 10*gamma/log(10)),R) + huber_pos((norm(x_var - anc(:,i)) - r(i)),R);
end
minimize obj
cvx_end

if (isnan(cvx_optval)) || (cvx_optval == +Inf) || (cvx_optval == -Inf)
    fail = true;
    return
end

x_est = x_var;

end

