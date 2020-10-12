function [alignments] = SPAL(N,m,sublen,matrix_profile,mp_index,...
            low_threshold,Mt,alignments)
    n_alignment = size(alignments,2);
    [ids, taken_subseqs] = get_taken_subseqs(N, sublen, alignments);
    RN = size(ids,1);
    if RN < 2
        return;
    end
    mplen = m-sublen+1;
    matching = zeros(RN,mplen);
    for i = 1:RN
        for j = 1:mplen
            for k = 1:n_alignment
                if taken_subseqs(ids(i),2*k-1) == 0
                    continue;
                end
                mxst = max(j,taken_subseqs(ids(i),2*k-1));
                mned = min(j+sublen-1,taken_subseqs(ids(i),2*k));
                count = mned-mxst+1;
                if count > 0
                    matching(i,j) = max(matching(i,j), count/sublen);
                end
            end
        end
    end
    candidates = cell(Mt,4);
    maxcor = zeros(Mt,1);
    mnidx = 1;
    c = 0;
    for i = 1:RN
        for j = 1:RN
            if i == j
                continue;
            end
            k = 1;
            while k <= mplen
                idx = mp_index(ids(i),ids(j),k)+1;
                if idx == 0
                    k = k+1;
                    continue;
                end
                cor = matrix_profile(ids(i),ids(j),k);
                if cor >= low_threshold
                    if matching(i,k) >= 0.5 || matching(j,idx) >= 0.5
                        k = k+1;
                        continue;
                    end
                    if cor >= maxcor(mnidx)
                        cand = {[i;j],k,idx,cor};
                        candidates(mnidx,:)= cand;
                        maxcor(mnidx) = cor;
                        if c+1 >= Mt
                            [mn,mnidx] = min(maxcor);
                            c = Mt;
                        else
                            c = c+1;
                            mnidx = mnidx+1;
                        end
                    end
                end
                k = k+1;
            end
        end
    end
    if c == 0
        clear candidates;
        return;
    end
    clear maxcor;
    candidates = sortrows(candidates(1:c,:),4,'descend');
    if RN > 2
        iteration = 2;
        Mt1 = floor(((Mt*(iteration+6))/(iteration+7))/(RN-iteration));
        c = min(Mt1, c);
        candidates = candidates(1:c,:);
        while iteration <= RN
            newcandidates = cell(Mt,4);
            cnew = 0;
            for i = 1:c
                cand = candidates(i,:);
                sz = size(cand{1},1);
                newset = zeros(sz+1,1);
                newset(1:sz) = cand{1};
                is_used = zeros(RN,1);
                for j = 1:sz
                    is_used(newset(j)) = 1;
                end
                k = newset(sz);
                idx = cand{3};
                startidx = cand{2};
                for j = 1:RN
                    if is_used(j) == 1
                        continue;
                    end
                    idx1 = mp_index(ids(k),ids(j),idx)+1;
                    if idx1 == 0
                        continue;
                    end
                    cor = matrix_profile(ids(k),ids(j),idx);
                    if cor >= low_threshold
                        if matching(j,idx1) >= 0.5
                            continue;
                        end
                        newset(sz+1) = j;
                        cnew = cnew+1;
                        newcandidates(cnew,:)= {newset,startidx,idx1,cand{4}+cor};
                    end
                end
            end
            if cnew > 0
                newcandidates = sortrows(newcandidates(1:cnew,:),4,'descend');
                iteration = iteration+1;
                if RN == iteration
                    candidates = newcandidates;
                    break;
                end
                Mt1 = floor(((Mt*(iteration+6))/(iteration+7))/(RN-iteration));
                c = cnew;
                c = min(Mt1, c);
                candidates = newcandidates(1:c,:);
            else
                break;
            end
        end
    end
    candidates = sortrows(candidates,4,'descend');
    mxidx = 1;
    n_alignment = n_alignment+1;
    alignments(n_alignment).sublen = sublen;
    alignments(n_alignment).lowt = low_threshold;
    cand = candidates(mxidx,:);
    sz = size(cand{1},1);
    tmp_ids = cand{1};
    idx = cand{2};
    act_ids = zeros(sz,1);
    subseqs = zeros(sz,1);
    act_ids(1) = ids(tmp_ids(1));
    subseqs(1) = idx;
    for i = 2:sz
        act_ids(i) = ids(tmp_ids(i));
        idx = mp_index(act_ids(i-1),act_ids(i),idx)+1;
        subseqs(i) = idx;
    end
    alignments(n_alignment).ids = act_ids;
    alignments(n_alignment).subseqs = subseqs;
    clear newcandidates;
    clear candidates;
end

function [ids,taken_subseqs] = get_taken_subseqs(N, cur_sublen, alignments)
    n_alignment = size(alignments,2);
    taken_subseqs = zeros(N,n_alignment*2);
    used = zeros(N,1);
    for i = 1:n_alignment
        sublen = alignments(i).sublen;
        ids = alignments(i).ids;
        subseqs = alignments(i).subseqs;
        n_ids = size(ids,1);
        for j = 1:n_ids
            taken_subseqs(ids(j),2*i-1) = subseqs(j);
            taken_subseqs(ids(j),2*i) = subseqs(j)+sublen-1;
            if sublen == cur_sublen
                used(ids(j)) = 1;
            end
        end
    end
    RN = N-sum(used);
    ids = zeros(RN,1);
    j = 0;
    for i = 1:N
        if used(i) == 0
            j = j+1;
            ids(j) = i;
        end
    end

end
