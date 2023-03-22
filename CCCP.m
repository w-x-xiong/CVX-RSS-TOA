function [x_est_vec, fail] = CCCP(d0, r, P0, P, gamma, anc, Nmax, epsilon)
%Concave-Convex Procedure (CCCP) for l1-Minimization
fail = false;
[H, L] = size(anc);

cnt = 0;

x_ini = zeros(H,1);

t_ini = zeros(2*L,1);

beta = zeros(L,1);

x_est_vec = [];

for i = 1:L
    beta(i) = 10^((P(i)-P0)/10/gamma);
end

x_est_vec = [x_est_vec,x_ini];

y_kth = [x_ini',t_ini']';

while cnt <= Nmax
    
    cvx_begin quiet
    
    variables y(H+2*L)
    
    expression obj_sum
    
    obj_sum = 0;
    
    for i = 1:2*L
        obj_sum = obj_sum + y(H+i);
    end
    
    minimize obj_sum

    subject to
    
    for i = 1:L
    
        (beta(i)/d0)*norm(y(1:H) - anc(:,i)) - 1 - (log(10)/(10*gamma))*y_kth(H+i) - (log(10)/(10*gamma))*[zeros(1,H),zeros(1,i-1),1,zeros(1,2*L-i)]*(y - y_kth) <= 0;
        
        1 - (log(10)/(10*gamma))*y(H+i) - (beta(i)/d0)*norm(y_kth(1:H) - anc(:,i)) - (beta(i)/d0)*[((y_kth(1:H)-anc(:,i))/norm(y_kth(1:H)-anc(:,i)))',zeros(1,2*L)]*(y - y_kth) <= 0;
        
        r(i) - y(H+L+i) - norm(y_kth(1:H) - anc(:,i)) - [((y_kth(1:H)-anc(:,i))/norm(y_kth(1:H)-anc(:,i)))',zeros(1,2*L)]*(y - y_kth) <= 0;
        
        norm(y(1:H) - anc(:,i)) - r(i) - y_kth(H+L+i) - [zeros(1,H),zeros(1,L+i-1),1,zeros(1,L-i)]*(y - y_kth) <= 0;
        
    end
    
    cvx_end
    
    if (isnan(cvx_optval)) || (cvx_optval == +Inf) || (cvx_optval == -Inf)
        fail = true;
        return
    end
    
    x_est_vec = [x_est_vec,y(1:H)];
    
    y_kth = y;
    
    cnt = cnt + 1;
    
    if (norm(x_est_vec(:,end) - x_est_vec(:,end-1))/(min(norm(x_est_vec(:,end)), norm(x_est_vec(:,end-1)))) <= epsilon)
        return
    end
    
end






end

