function [kAlignments, n_alignment] = MultiPal(N, m, K, L, Mt)
    alignments = struct('sublen',{},'lowt',{},'ids',{},'subseqs',{});
    n_alignment = 0;
    S = length(L);
    for i = 1:S
        sublen = L(i);
        mplen = m-sublen+1;
        matrix_profile = zeros(N, N, mplen);
        mp_index = zeros(N, N, mplen);
        for j = 1:N
            for k = 1:N
                if j ~= k
                    filename = strcat(int2str(sublen),'/mp_dist',int2str(j),'_',int2str(k),'.txt');
                    if isfile(filename)
                        dist = importdata(filename);
                        dist = 1-((dist.^2)./(2*sublen));
                        if length(dist) < mplen
                            for l = length(dist)+1:mplen
                                dist(l) = NaN;
                            end
                        end
                        matrix_profile(j,k,:) = dist;
                        filename = strcat(int2str(sublen),'/mp_index',int2str(j),'_',int2str(k),'.txt');
                        data = importdata(filename);
                        mp_index(j,k,:) = data;
                    else
                        fprintf('File %s not found\n',filename);
                    end
                end
            end
        end
        [low_threshold] = get_threshold(matrix_profile,90);
        low_threshold = round(low_threshold,3);
        %fprintf('Running for sublen %d\n', sublen);
        while 1
            alignments = SPAL(N,m,sublen,matrix_profile,mp_index,...
                low_threshold,Mt,alignments);
            if n_alignment == length(alignments)
                break;
            end
            n_alignment = n_alignment+1;
        end
        clear matrix_profile;
        clear mp_index;
    end
    [kAlignments, n_alignment] = combine_alignments(N, alignments, K);
end
