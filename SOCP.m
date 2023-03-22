function [x_est, fail] = SOCP(d0, r, P0, P, gamma, anc)
%Directly relaxing the epigraph-form problem to an SOC program

fail = false;

[H, L] = size(anc);
x_est = zeros(H,1);
beta = zeros(L,1);

for i = 1:L
    beta(i) = 10^((P(i)-P0)/10/gamma);
end

    cvx_begin quiet
    
    variables y(H+2*L) d(L)
    
    expression obj_sum
    
    obj_sum = 0;
    
    for i = 1:2*L
        obj_sum = obj_sum + y(H+i);
    end
    
    minimize obj_sum

    subject to
    
    for i = 1:L
    
        abs(10*gamma*(beta(i)*d(i)/d0 - 1)/log(10)) <= y(H+i);
        
        abs(r(i) - d(i)) <= y(H+i+L);
        
        norm(y(1:H) - anc(:,i)) <= d(i);
        
    end
    
    cvx_end
    
    if (isnan(cvx_optval)) || (cvx_optval == +Inf) || (cvx_optval == -Inf)
        fail = true;
        return
    end
    
    x_est = y(1:H);

end

