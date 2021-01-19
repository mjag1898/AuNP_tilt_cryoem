% Kabsch algorithm for rotating mobile to target to minimize rmsd
% Assumes target and mobile have been centered at the origin

function [R, rmsd, rmsds] = kabsch(target, mobile, num_points)
    [nmols, ~] = size(target);
    Rs = [];
    rmsds = [];
    if nargin < 3
        num_points = nmols;
    end
    
    selections = combnk(1:nmols, num_points);
    num_selections = nchoosek(nmols,num_points);
    for i = 1:num_selections
        A = mobile(selections(i,:),:)'*target(selections(i,:),:);
        R = ((A'*A)^0.5)*inv(A);
        rmsd = norm((R*mobile(selections(i,:),:)') - target(selections(i,:),:)', 'fro')/sqrt(num_points);
        rmsds = [rmsds rmsd];
        Rs = cat(3, Rs, R);
    end
    [~, ind] = min(rmsds);
    R = Rs(:,:,ind);
    rmsd = rmsds(ind);
end