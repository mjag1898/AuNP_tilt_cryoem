% After finalizing the individual 3D structures, do a final pass of
% aligning all structures to each other.

% Use model coordinates as center for simulated data or use a real
% structure as center for real data
if real == 1
    center = structures(:,:,best_mean);
elseif real == 0
    center = [[ 12.234589,  -6.453163, -98.69324 ]; ...
       [  4.263031, -31.591179,  71.67932 ]; ...
       [ 21.443192,  12.404953,   6.341053]; ...
       [-37.940796,  25.639404,  20.672852]];
end

fit_structures = zeros(size(structures));
[~, ~, n] = size(structures);
configs = perms(1:4);

for iter = 1:8
    for i = 1:n
        best_rmsd = 100000;
        best_R = 0;
        best_config = 0;
        for k = 1:24
            [R, rmsd, ~] = kabsch(center, structures(configs(k,:), :, i));
            if rmsd < best_rmsd
                best_rmsd = rmsd;
                best_R = R;
                best_config = configs(k,:);
            end
        end

        fit_structures(:,:,i) = structures(best_config, :, i)*best_R';
    end
    center = mean(fit_structures, 3);
end

diffs = fit_structures - center;
structures = fit_structures;


vars = {'best_config','best_R','best_rmsd','configs','i','k','iter'...
    'n','R','rmsd','vars','center','diffs','fit_structures'};
clear(vars{:})