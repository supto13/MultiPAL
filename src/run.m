clear;
load('../data/data.mat');
N = 5;
m = 3000;
k = 3;
t_k = k;
low_b = 5;
up_b = 15;
iteration = 1;
while iteration <= 10
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
    [kAlignments, n_alignment] = MultiPal(N, m, k, L, 500000);
    if n_alignment == k
        save(strcat(int2str(k),'-alignments.mat'),'kAlignments','n_alignment');
        for i = 1:n_alignment
            plot_patterns(kAlignments(i), i, data, 1);
            pause;
            close all;
            display_alignments(kAlignments(i), i, data, 1);
            pause;
            close all;
        end
        break;
    end
    iteration = iteration+1;
    low_b = low_b-1;
    up_b = up_b+1;
    t_k = t_k+2;
end
if iteration > 5
    fprintf('No %d-alignment found\n', k);
end

