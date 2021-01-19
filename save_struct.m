% Save set of reconstructed 3D structures

diffs = structures - mean(structures, 3);
rmsds = permute(sqrt(mean(sum((diffs .^ 2), 2), 1)), [3 2 1]);
[~, I] = sort(rmsds);
structures = structures(:, :, I);
center = mean(structures, 3);
%center = [[ 12.234589,  -6.453163, -98.69324 ]; ...
%       [  4.263031, -31.591179,  71.67932 ]; ...
%       [ 21.443192,  12.404953,   6.341053]; ...
%       [-37.940796,  25.639404,  20.672852]];
structures = cat(3, center, structures);

structures_struct.structs = structures;
save('real_structs.mat','-struct','structures_struct')


vars = {'diffs','rmsds','I','structures','disagreement_mat','tree','v',...
    'i','j','structures_struct','vars','rmsd'...
    };
clear(vars{:})
