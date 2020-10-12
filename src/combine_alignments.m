function [comb_alignments, n_alignment] = combine_alignments(N, alignments,K)
    comb_alignments = struct('sublens',{},'lowts',{},'ids',{},'subseqs',{});
    n_alignment = size(alignments,2);
    for i = 1:K
        C = combnk(1:n_alignment,i);
        len = size(C,1);
        mxcount = 0;
        mxsublens = zeros(i,1);
        mxlow = zeros(i,1);
        mxids = zeros(N,1);
        mxsubseqs = zeros(N,i);
        mn_mutual_info = 1000000000;
        for j = 1:len
            alignment_ids = C(j,:);
            visited = zeros(N,1);
            all_subseqs = zeros(N,i);
            tmp_sublens = zeros(i,1);
            tmp_low = zeros(i,1);
            for k = 1:i
                aid = alignment_ids(k);
                ids = alignments(aid).ids;
                tmp_subseqs = alignments(aid).subseqs;
                tmp_sublens(k) = alignments(aid).sublen;
                tmp_low(k) = alignments(aid).lowt;
                n_ids = size(ids,1);
                for l = 1:n_ids
                    visited(ids(l)) = visited(ids(l))+1;
                    all_subseqs(ids(l),k) = tmp_subseqs(l);
                end
            end
            cnt = 0;
            tmp_ids = zeros(N,1);
            tmp_subseqs = zeros(N,i);
            for l = 1:N
                if visited(l) == i
                    cnt = cnt+1;
                    tmp_ids(cnt) = l;
                    tmp_subseqs(cnt,:) = all_subseqs(l,:);
                end
            end
            tmp_ids = tmp_ids(1:cnt);
            tmp_subseqs = tmp_subseqs(1:cnt,:);
            if cnt >= 2
                if cnt > mxcount
                    mxcount = cnt;
                    mxsublens = tmp_sublens;
                    mxlow = tmp_low;
                    mxids = tmp_ids;
                    mxsubseqs = tmp_subseqs;
                    mn_mutual_info = compute_mutual_info(mxcount,i,mxsublens,mxsubseqs);
                elseif cnt == mxcount
                    mutual_info = compute_mutual_info(cnt,i,tmp_sublens,tmp_subseqs);
                    if mutual_info < mn_mutual_info || (mutual_info == mn_mutual_info &&...
                            max(tmp_sublens) > max(mxsublens))
                        mxsublens = tmp_sublens;
                        mxlow = tmp_low;
                        mxids = tmp_ids;
                        mxsubseqs = tmp_subseqs;
                        mn_mutual_info = mutual_info;
                    end
                end
            end
        end
        
        if mxcount > 0
            comb_alignments(i).sublens = mxsublens;
            comb_alignments(i).lowts = mxlow;
            comb_alignments(i).ids = mxids;
            comb_alignments(i).subseqs = mxsubseqs;
        else
            n_alignment = i-1;
            return;
        end
    end
    n_alignment = K;
end
function [mutual_info] = compute_mutual_info(n_id, n_subseq, sublens, subseqs)
    mutual_info = 0;
    cnt = 0;
    for i = 1:n_id
        for j = 1:n_subseq
            for k = j+1:n_subseq
                st = subseqs(i,j);
                ed = st+sublens(j)-1;
                st1 = subseqs(i,k);
                ed1 = st1+sublens(k)-1;
                stmx = max(st,st1);
                mned = min(ed,ed1);
                if stmx <= mned
                    mutual_info = mutual_info+(2*(mned-stmx+1))/(sublens(j)+sublens(k));
                end
                cnt = cnt+1;
            end
        end
    end
    if cnt > 0
        mutual_info = mutual_info/cnt;
    end
end
