% Align all structures to each other
% Start with a single 3D structure and align structures to the mean one at a time
% starting with those closest to the mean, updating the mean after each structure. 
% Test all possible starting structures and mark the best (lowest total RMSD to mean)

[~, ~, n] = size(structures);
mean_structures = zeros(4, 3, n);
configs = perms(1:4);

for i = 1:n
    curr_mean = structures(:,:,i);
    rmsds = 10000*ones(n,1);
    rmsds(i) = 0;
    
    for j = setdiff(1:n, i)
        to_add = structures(:,:,j);
        for k = 1:24
            [R, rmsd, ~] = kabsch(curr_mean, to_add(configs(k,:), :));
            if rmsd < rmsds(j)
                rmsds(j) = rmsd;
            end
        end
    end
    
    [~, add_order] = sort(rmsds);
    structures_sort = structures(:,:,add_order);
    
    for j = 2:n
        to_add = structures_sort(:,:,j);
        best_rmsd = 10000;
        best_config = 0;
        best_R = 0;
        for k = 1:24
            [R, rmsd, ~] = kabsch(curr_mean, to_add(configs(k,:), :));
            if rmsd < best_rmsd
                best_rmsd = rmsd;
                best_config = configs(k,:);
                best_R = R;
            end
        end
        
        to_add = to_add(best_config, :);
        to_add = to_add * best_R';
        curr_mean = (1/j)*((j-1)*curr_mean + to_add);
    end
    
    mean_structures(:,:,i) = curr_mean;
end

diffs = sum((mean_structures - structures).^2, 1);
diffs = sum(diffs, 2);
diffs = sqrt(diffs);
diffs = diffs(:);
[~, best_mean] = min(diffs);

vars = {'best_perm','best_R','best_rmsd','configs','curr_mean','i','j',...
    'k','n','N','R','rmsd','to_add','vars','add_order',...
    'best_config','rmsds','structures_sort','diffs','mean_structures'};
clear(vars{:})
