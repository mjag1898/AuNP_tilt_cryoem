% Set to 1 for real data, 0 for sim data. Specific to our dataset and file setup
real = 1;
% Threshold for least-squares objective
rmsd_thres = 5;
% Pixel size in Angstroms
pixel_size = 1.285;


if real == 1
    % Settings for where to find images
    file_nums = [1:49,79:123];
    num_mols = length(file_nums);
    tifFileformat = {'part0deg-%03d.tif','part-45deg-%03d.tif'};
    image_dir = 'EM_images_real';
    % Settings for off-axis angle
    l = -cos(86.5*pi/180);
    m = cos(3.5*pi/180);
elseif real == 0
    % Settings for where to find images
    file_nums = 1:141;
    num_mols = length(file_nums);
    tifFileformat = {'stack_0deg-%03d.tif', 'stack_neg45deg-%03d.tif'};
    image_dir = 'EM_images_sim';
    % Settings for off-axis angle
    l = -cos(90*pi/180);
    m = cos(0*pi/180);
end

% Load all images and find circle centers in each one
% Tag images with less than 4 circles for removal 
% Tag images with more than 4 circles for further analysis
% See code for function find_Au
allcenters=[];
allRcenters=[];
for k = 1:num_mols
    tifFilename = sprintf(char(tifFileformat(1)), file_nums(k));
    fullFileName = fullfile(image_dir, tifFilename);
    
    imageFileName = sprintf('stack_0deg_P-%d.tif', k);
    centers = find_Au(fullFileName, 0, imageFileName);
    
    [rows,~] = size(centers);
    % Scale
    centers(:, 1:3) = pixel_size * centers(:, 1:3);
    centers = [k*ones(rows,1) centers];
    allcenters=[allcenters; centers];
    
    if k ==1
        w = warning('query','last');
        id = w.identifier;
        warning('off',id);
    end

    tifFilename = sprintf(char(tifFileformat(2)), file_nums(k));
    fullFileName = fullfile(image_dir, tifFilename);
    
    imageFileName = sprintf('stack_neg45deg_P-%d.tif', k);
    centers = find_Au(fullFileName, 0, imageFileName);
    
    [rows,~] = size(centers);
    % Scale
    centers(:, 1:3) = pixel_size * centers(:, 1:3);
    centers = [k*ones(rows,1) centers];
    allRcenters=[allRcenters; centers];       
end

% Remove all tilt pairs where either image has less than 4 circles
allcenters = allcenters(allcenters(:,5) ~= 2,:);
allRcenters = allRcenters(allRcenters(:,5) ~= 2,:);
allcenters = allcenters(ismember(allcenters(:,1), allRcenters(:,1))...
    | allcenters(:,1) == 0,:);
allRcenters = allRcenters(ismember(allRcenters(:,1), allcenters(:,1)),:);

% Make arrays for circles in untilted and tilted images the same shape
for k = 1:num_mols
    orig = allcenters(allcenters(:,1) == k,:);
    rotate = allRcenters(allRcenters(:,1) == k,:);
    [orig_rows,~] = size(orig);
    [rotate_rows,~] = size(rotate);

    num_rows = max(orig_rows, rotate_rows);
    new_orig = nan*ones(num_rows, 5);
    new_rotate = nan*ones(num_rows, 5);
    new_orig(:, 1) = k;
    new_rotate(:, 1) = k;
    new_orig(1:orig_rows, :) = orig;
    new_rotate(1:rotate_rows, :) = rotate;

    allcenters = [allcenters(allcenters(:,1) < k,:); new_orig; allcenters(allcenters(:,1) > k,:)];
    allRcenters = [allRcenters(allRcenters(:,1) < k,:); new_rotate; allRcenters(allRcenters(:,1) > k,:)];
end

CombM=[allcenters allRcenters];
CombM = CombM(:,[1:3,7:8]);

% Construct coordinate rotation matrix between untilted and tilted image
n = 0;
deg = 45;
s=deg*pi/180;
R11 = l*l*(1-cos(s))+cos(s);
R12 = (m*l*(1-cos(s))-n*sin(s));
R13 = (n*l*(1-cos(s))+m*sin(s));
R21 = (l*m*(1-cos(s))+n*sin(s));
R22 = (m*m*(1-cos(s))+cos(s));
R23 = (n*m*(1-cos(s))-l*sin(s));
R31 = (l*n*(1-cos(s))-m*sin(s));
R32 = (m*n*(1-cos(s))+l*sin(s));
R33 = (n*n*(1-cos(s))+cos(s));

Rmatrix=[R11,R12, R13;R21, R22, R23;R31, R32, R33];  
drop_mat = [1 0 0; 0 1 0];
num_left = length(unique(allcenters(:, 1)));
Solution_Mat = zeros(4*num_left, 3);
Solution_Scores = zeros(num_left,1);

% Matrix for use in least squares calculation
B = [drop_mat zeros(2, 9); drop_mat*Rmatrix zeros(2,9); ...
    zeros(2,3) drop_mat zeros(2,6); zeros(2,3) drop_mat*Rmatrix zeros(2,6); ...
    zeros(2,6) drop_mat zeros(2,3); zeros(2,6) drop_mat*Rmatrix zeros(2,3); ...
    zeros(2,9) drop_mat; zeros(2,9) drop_mat*Rmatrix];

more_than_4 = zeros(num_left, 1);
mul_good_corr = zeros(num_left, 1);

% 1) For each tilt-pair, iterate through every selection of 4 particles 
% and every correspondence between particles in the two images
% 2) For each selection and correspondence, compute y-dev and the 3D reconstruction
% along with its least-squares objective value
% 3) Store best 3D structure for each tilt-pair and also mark which ones had multiple
% good solutions and require further analysis.
count = 1;
check = zeros(num_left, 1);
check_mols = [];
all_y_dev = {};
all_errors = {};
for mol = 1:num_mols
    orig = CombM(CombM(:, 1) == mol, 2:3);
    rotate = CombM(CombM(:, 1) == mol, 4:5);
    n_orig = sum(~isnan(orig(:, 1)));
    n_rotate = sum(~isnan(rotate(:, 1)));
    
    if (n_orig > 4 || n_rotate > 4)
        more_than_4(count) = 1;
    end
    
    if (n_orig == 0 || n_rotate == 0)
        continue;
    end
    
    orig_comb = nchoosek(1:n_orig,4);
    rot_perms = perms(1:n_rotate);
    rot_perms = rot_perms(:, 1:4);
    [num_comb, ~] = size(orig_comb);
    [num_perms, ~] = size(rot_perms);
    
    errors = zeros(num_comb*num_perms, 1);
    y_dev = zeros(num_comb*num_perms, 1);
    comb_list = zeros(num_comb*num_perms, 1);
    perm_list = zeros(num_comb*num_perms, 1);
    
    for comb = 1:num_comb
        test_orig = orig(orig_comb(comb, :), :);
        test_orig = test_orig - mean(test_orig, 1);
        
        for perm = 1:num_perms
            test_rotate = rotate(rot_perms(perm, :), :);
            test_rotate = test_rotate - mean(test_rotate, 1);
            y_dev(num_perms*(comb-1)+perm) = ...
                sqrt(sum((test_rotate(:,2)-test_orig(:,2)).^2));
            
            c = reshape([test_orig, test_rotate]', [16, 1]);
            sol = B\c;
            errors(num_perms*(comb-1)+perm) = norm(B*sol - c) / sqrt(8);
            
            comb_list(num_perms*(comb-1)+perm) = comb;
            perm_list(num_perms*(comb-1)+perm) = perm;
        end
    end
    all_y_dev = [all_y_dev y_dev];
    all_errors = [all_errors errors];
    if (min(errors) < rmsd_thres)
        [~, I] = sort(y_dev);
        errors = errors(I);
        comb_list = comb_list(I);
        perm_list = perm_list(I);
        I = find(errors < rmsd_thres, 1);
        if sum(errors < rmsd_thres & sort(y_dev) < 6) > 1
            mul_good_corr(count) = 1;
            check(count) = 1;
            check_mols = [check_mols; mol];
        end
    else
        [~, I] = min(errors);
    end
    comb = comb_list(I);
    perm = perm_list(I);
    best_orig = orig(orig_comb(comb, :), :);
    best_orig = best_orig - mean(best_orig, 1);
    best_rotate = rotate(rot_perms(perm, :), :);
    best_rotate = best_rotate - mean(best_rotate, 1);

    c = reshape([best_orig, best_rotate]', [16, 1]);
    best_sol = reshape(B\c, [3 4])';
    
    Solution_Scores(count) = errors(I);
    Solution_Mat(4*count - 3:4*count, :) = best_sol;
    count = count + 1;
end

structures = reshape(Solution_Mat, 4, [], 3);
structures = permute(structures, [1 3 2]);
structures = structures(:, :, Solution_Scores < rmsd_thres);
check = check(Solution_Scores < rmsd_thres);
check = find(check);

warning('on',id);
vars = {'vars','allcenters','allRcenters','c','centers','comb','count',...
    'deg','drop_mat','fullFileName','k','l','m','n','new_orig','new_rotate',...
    'orig','orig_rows','R11','R12','R13','R21','R22','R23','R31','R32','R33',...
    'rotate','rotate_rows','row','rows','s','sol','mol','n_orig',...
    'tifFilename','num_mols','tifFileformat','image_dir','w','id',...
    'rot','rot_perms','test_rotate','n_rotate','y_dev','errors','comb_list',...
    'perm_list','num_comb','num_left','num_perms','num_rows','orig_comb',...
    'perm','test_orig','best_orig','best_rotate','best_sol','I',...
    'Solution_Mat','all_errors','all_y_dev','keep','file_nums',...
    'imageFileName'};
    
clear(vars{:})

%{
hold on;
histogram(sol_scores_real, [0:0.418:8.8])
histogram(sol_scores_sim, [0:0.418:8.8])
legend({'Real data', 'Simulated data'}, 'FontSize', 14)
set(gca, 'FontSize', 14)
Ang = char(197);
xlabel(strcat('Coordinate to image RMSD (', Ang, ')'))
title('Histogram of 3D coordinate quality')
ylabel('Count')
%}
