load('accelY.mat');
for i = 1:1
    file_dir = 'output_merge/pattern';
    filename = strcat(file_dir,int2str(i),'/1.txt');
    fprintf('%s\n', filename);
    if isfile(filename)
        file_dir = strcat(file_dir,int2str(i),'/');
        XX = zeros(20,1);
        for j = 1:20
            filename = strcat(file_dir,int2str(j),'.txt');
            if isfile(filename)
                plot_patterns(file_dir, j, accelY, 1);
                pause;
                close all;
                fprintf('Combination %d\n', j);
                display_alignments(file_dir, j, accelY, 1);
                pause;
                close all;
            else
                break;
            end
        end
    else
        break;
    end
end
