% Make 3D plot of reconstructed 3D structures showing radius of AuNP

function plot_mols(data, sz, alpha, edge)
    if nargin < 2
        sz = 8;
    end
    if nargin < 3
        alpha = 1;
    end
    if nargin < 4
        edge = 0;
    end
    [~, ~, sze] = size(data);
    hold on;
    colors = [0 1 0; 0 0 0; 0 0 1; 1 0 0];
    for i = [1:sze]
        if (i == 1)
            color = [0 0 0];
        else
            color = hsv2rgb([i/sze, 1, 1]);
        end
        x = data(1, :, i);
        y = data(2, :, i);
        z = data(3, :, i);
        disp = 0.2;
        if edge == 1
            scatter3(x, y, z, sz, colors, 'filled');
            %text(x+disp, y+disp, z+disp, cellstr(num2str(i)),'Fontsize',7);
        elseif edge == 0
            for j = 1:4
                [xs, ys, zs] = sphere;
                hSurface = surf(sz*xs+x(j), sz*ys+y(j), sz*zs+z(j));
                set(hSurface, 'Facecolor',colors(j, :),'FaceAlpha',alpha,...,
                    'FaceLighting','gouraud','EdgeColor','none');
            end
        end
    end
    rotate3d on;
    axis equal;
    grid on;
    snapnow;
end
