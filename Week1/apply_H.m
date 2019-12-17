function [transformedImage] = apply_H(image, H)
% Apply a homography to an image
%


image = double(image);  
s = size(image);


boundaries = [1, s(2), s(2), 1; 1, 1, s(1), s(1),; 1, 1, 1, 1];
tl=H*boundaries(:,1);
tl=tl/tl(3);
tr=H*boundaries(:,2);
tr=tr/tr(3);
br=H*boundaries(:,3);
br=br/br(3);
bl=H*boundaries(:,4);
bl=bl/bl(3);
point_boundaries = [tl';tr';br';bl'];
new_boundaries = ceil([max(point_boundaries(:,1)), max(point_boundaries(:,2))]);
transformedImage = uint8(zeros(new_boundaries(2), new_boundaries(1), 3));

for x=1: new_boundaries(1)
    for y=1: new_boundaries(2)
        xy = [x, y, 1]';
        original_point = H\xy;
        original_point = round(original_point/original_point(3));
        
        if (original_point(1) < s(2) && original_point(2) < s(1))
            if (original_point(1) > 0 && original_point(2) > 0)
                transformedImage(y, x, :) = image(original_point(2), original_point(1),:);
            end
        end
    end
end


