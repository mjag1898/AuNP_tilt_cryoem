% Find circle centers in a single tif image

function centers = find_Au(fullFileName, noise, imageFileName)
    if (exist(fullFileName, 'file'))
        imageData = imread( fullFileName );
    else
        FinishMessage = sprintf('Warning: image file does not exist:\n%s', fullFileName);
        uiwait(warndlg(FinishMessage));
    end
    
    Igray=mat2gray(imageData);
    RevOriI=1-Igray;
    WienerI=wiener2(Igray,[5,5]);
    AdjustI=imadjust(WienerI);
    RevI=1-AdjustI;
    level=0.95;
    Ithresh=im2bw(RevI,level);
    [circle_centers, radii]=imfindcircles(Ithresh, [5,10], 'Sensitivity',0.95);
    circle_centers = circle_centers + noise * randn(size(circle_centers));
    
    % Write figure file to show how AuNPs are recognized in EM images
    %imshow(RevOriI);
    %h=viscircles(circle_centers,radii);
    %saveas(gcf,imageFileName);
    
    centers_radii=[circle_centers radii];
    
    sorted_centers=sortrows(centers_radii, 2);
    [rows, ~] = size(sorted_centers);

    % If 4 or less found, center the coordinates at the origin.
    % If more than 4, do not do this because we have to test different choice of 4 later
    if (rows < 5)
        col = 1;
        SumX = 0;
        for row = 1 : rows
            SumX = SumX + sorted_centers(row, col);
        end
        CenterOfX=SumX/rows;

        col = 2;
        SumY = 0;
        for row = 1 : rows
            SumY = SumY + sorted_centers(row, col);
        end
        CenterOfY=SumY/rows;
         
        moved_centers=sorted_centers;
        col = 1;
        for row = 1 : rows
            moved_centers(row, col) = moved_centers(row, col) - CenterOfX;
        end
        col = 2;
        for row = 1 : rows
            moved_centers(row, col) = moved_centers(row, col) - CenterOfY;
        end
    else
        moved_centers = sorted_centers;
    end
   
    % If less than 4 found, flag for removal
    % If more than 4 found, flag for further analysis

    if rows > 4
        remove = ones(rows,1);
    elseif rows < 4
        remove = 2*ones(rows,1);
    else
        remove = zeros(rows,1);
    end
    centers = [moved_centers remove];
end

