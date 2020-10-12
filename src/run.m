clear;
N = 5;
m = 3000;
for K = 2:10
    t_k = K;
    low_b = 5;
    up_b = 15;
    iteration = 1;
    while iteration <= 5
        L = zeros(t_k,1);
        interval = (up_b-low_b)/(t_k-1);
        for i = 1:t_k
            L(i) = ceil(((low_b+interval*(i-1))/100)*m);
        end
        tmpL = zeros(t_k,1);
        for i = 1:t_k
            tmpL(i) = L(t_k-i+1);
        end
        L = tmpL;
        [kAlignments, n_alignment] = MultiPal(N, m, K, L, 500000);
        if n_alignment == K
            fprintf('%d %d\n', K, t_k);
            break;
        end
        iteration = iteration+1;
        low_b = low_b-1;
        up_b = up_b+1;
        t_k = t_k+2;
    end
    if iteration > 5
        fprintf('%d %d\n', K, t_k-2);
    end
end
