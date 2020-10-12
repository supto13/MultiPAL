function [] = plot_patterns(alignment, n_comb, data, is_normalized)
    sublens = alignment.sublens;
    low_thresholds = alignment.lowts;
    up_thresholds = alignment.upts;
    ids = alignment.ids;
    subseqs = alignment.subseqs;
    n_ids = size(ids,1);
    colors = ["#006400";"#CD853F";"#FF6EEC";"#0000CD";"#008B8B";"#FF0000";"#DAA520"];
    M = size(data, 2);
    XX = 1:M;
    legend_labels = strings(n_comb,1);
    lgd = zeros(1, n_comb);
    if is_normalized == 1
        for i = 1:n_ids
            data(ids(i),:) = zNorm(data(ids(i),:));
        end
    end
    shift = 0;
    for i = 1:n_comb
        mx = -1000000000;
        id = ids(1);
        subid = subseqs(1,i);
        if i ~= 1
            shift = shift+abs(min(0,min(data(id,subid:subid+sublens(i)-1))));
        end
        legend_label = strcat('sublen=',int2str(sublens(i)),',low thrshld=',...
            num2str(low_thresholds(i)),',up thrsld=',num2str(up_thresholds(i)));
        legend_labels(i) = legend_label;
        lgd(i) = plot(XX(1:sublens(i)),shift+data(id,subid:subid+sublens(i)-1),...
            'Color',colors(i));
        hold on;
        mx1 = max(data(id,subid:subid+sublens(i)-1));
        mn1 = min(0,min(data(id,subid:subid+sublens(i)-1)));
        mx = max(mx, mx1-mn1);
        shift = shift+mx+1.0;
    end
    hold off;
    figname = strcat('# of aligned timeseries=',int2str(n_ids),',# of combination=',...
        int2str(n_comb));
    title(figname);
    legend(lgd, cellstr(legend_labels));
end