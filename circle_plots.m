% Generate 3D plots of recontstructions 

hold on;
if real == 1
    mean_struct = mean(structures, 3);
elseif real == 0      
    mean_struct = [[ 12.234589,  -6.453163, -98.69324 ]; ...
           [  4.263031, -31.591179,  71.67932 ]; ...
           [ 21.443192,  12.404953,   6.341053]; ...
           [-37.940796,  25.639404,  20.672852]];
end
plot_mols(permute(structures / 10, [2 1 3]), 8, 1, 1);
plot_mols(mean_struct' / 10, 0.8, 0.2, 0);
set(gca, 'FontSize', 20)
xlabel('X (nm)')
ylabel('Y (nm)')
zlabel('Z (nm)')