function [] = display_alignments(alignment, n_comb, data, is_normalized)
    sublens = alignment.sublens;
    low_thresholds = alignment.lowts;
    up_thresholds = alignment.upts;
    ids = alignment.ids;
    subseqs = alignment.subseqs;
    n_ids = size(ids,1);
    shift = 0;
    colors = ["#ADFF2F";"#FFE812";"#FF00FF";"#87CEFA";"#00FFFF";"#FF8383";"#FFFF00"];
    darkcolors = ["#006400";"#CD853F";"#8B008B";"#0000CD";"#008B8B";"#FF0000";"#DAA520"];
    M = size(data,2);
    XX = 1:M;
    legend_labels = strings(n_comb,1);
    lgd = zeros(1, n_comb);
    for i = 1:n_ids
        if i > 10
            break;
        end
        
        id = ids(i);
        if is_normalized == 1
            data(id,:) = zNorm(data(id,:));
        end
        if i ~= 1
            shift = shift+abs(min(0,min(data(id,:))));
        end
        mx = abs(max(data(id,:)));
        data(id,:) = shift+data(id,:);
        plot(XX,data(id,:),'Color','#696969');
        hold on;
        for j = 1:n_comb
            subid = subseqs(i,j);
            if i == 1
                legend_label = strcat('sublen=',int2str(sublens(j)),',low threshld=',...
                num2str(low_thresholds(j)),',up threshld=',num2str(up_thresholds(j)));
                legend_labels(j) = legend_label;
                lgd(j) = plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',colors(j),'LineWidth',8);
                plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',darkcolors(j),'LineWidth',2);
            else
                plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',colors(j),'LineWidth',8);
                 plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',darkcolors(j),'LineWidth',2);
            end
            hold on;
        end
        shift = shift+mx+1.0;
    end
    hold off;
    xlabel('Timeline');
    if is_normalized == 1
        ylabel('Z-Normalized value (Shifted)');
    else
        ylabel('Actual value (Shifted)');
    end
    
    figname = strcat('# of aligned timeseries=',int2str(n_ids),',# of combination=',...
        int2str(n_comb));
    title(figname);
    legend(lgd, cellstr(legend_labels));
end
