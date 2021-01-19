% Picks the 3D structure for particles that had multiple high quality options
% Chooses the structure that matches the mean of all structures the best by RMSD.
% Use after initial pass of aligning all structures to each other.

mean_struct = mean(structures, 3);
for num = 1:length(check)
    mol = check_mols(num);
    orig = CombM(CombM(:, 1) == mol, 2:3);
    rotate = CombM(CombM(:, 1) == mol, 4:5);
    n_orig = sum(~isnan(orig(:, 1)));
    n_rotate = sum(~isnan(rotate(:, 1)));
    if (n_orig == 0 || n_rotate == 0)
        continue;
    end
    
    orig_comb = nchoosek(1:n_orig,4);
    orig_perm = perms(1:4);
    rot_perms = perms(1:n_rotate);
    rot_perms = rot_perms(:, 1:4);
    [num_comb, ~] = size(orig_comb);
    [num_perms, ~] = size(rot_perms);
    
    errors = zeros(num_comb*num_perms*24, 1);
    rmsds = zeros(num_comb*num_perms*24, 1);
    comb_list = zeros(num_comb*num_perms*24, 1);
    perm_orig_list = zeros(num_comb*num_perms*24, 1);
    perm_list = zeros(num_comb*num_perms*24, 1);
    
    for comb = 1:num_comb
        test_orig_comb = orig(orig_comb(comb, :), :);
        test_orig_comb = test_orig_comb - mean(test_orig_comb, 1);
        
        for perm_1 = 1:24
            test_orig = test_orig_comb(orig_perm(perm_1, :), :);
            
            for perm = 1:num_perms
                index = (comb-1)*24*num_perms + (perm_1-1)*num_perms + perm;
                test_rotate = rotate(rot_perms(perm, :), :);
                test_rotate = test_rotate - mean(test_rotate, 1);

                c = reshape([test_orig, test_rotate]', [16, 1]);
                sol = B\c;
                errors(index) = norm(B*sol - c) / sqrt(8);
                sol = reshape(sol, [3 4])';
                [~, rmsd, ~] = kabsch(mean_struct, sol);
                rmsds(index) = rmsd;

                comb_list(index) = comb;
                perm_orig_list(index) = perm_1;
                perm_list(index) = perm;
            end
        end
    end
    rmsds(errors > rmsd_thres) = Inf;
    [~, I] = min(rmsds);
    comb = comb_list(I);
    perm = perm_list(I);
    perm_orig = perm_orig_list(I);
    best_orig = orig(orig_comb(comb, :), :);
    best_orig = best_orig - mean(best_orig, 1);
    best_orig = best_orig(orig_perm(perm_orig, :),:);
    best_rotate = rotate(rot_perms(perm, :), :);
    best_rotate = best_rotate - mean(best_rotate, 1);
    c = reshape([best_orig, best_rotate]', [16, 1]);
    best_sol = reshape(B\c, [3 4])';
    [R, ~, ~] = kabsch(mean_struct, best_sol);
    structures(:, :, check(num)) = best_sol * R';
end

vars = {'B','best_orig','best_rotate','best_sol','c','check',...
    'check_mols','comb','comb_list','CombM','errors','I','index','mean_struct',...
    'mol','n_orig','n_rotate','num','num_comb','num_perms','orig','orig_comb',...
    'orig_perm','perm','perm_1','perm_list','perm_orig','perm_orig_list',...
    'rmsds','rot_perms','rotate','sol','test_orig','test_orig_comb','test_rotate',...
    'vars','rmsd','R','rmsd_thres'};
clear(vars{:})