function [transformedImage, x_min, y_min, x_max, y_max] = apply_H(image, H)
% Apply a homography to an image
%

    image = double(image);  
    s = size(image);

    topleft = H * [0; 0; 1];        
    topright = H * [0; s(2); 1];     
    bottomleft = H * [s(1); 0; 1];     
    bottomright = H * [s(1); s(2); 1];  

    topleft = topleft / topleft(3);
    topright = topright / topright(3);
    bottomleft = bottomleft / bottomleft(3);
    bottomright = bottomright / bottomright(3);

    x_min = round(min([topleft(1), topright(1), bottomleft(1), bottomright(1)]));
    y_min = round(min([topleft(2), topright(2), bottomleft(2), bottomright(2)]));
    x_max = round(max([topleft(1), topright(1), bottomleft(1), bottomright(1)]));
    y_max = round(max([topleft(2), topright(2), bottomleft(2), bottomright(2)]));
    out_size = [abs(x_max-x_min), abs(y_max-y_min), 3];
    %out_size = s;

    x = zeros(out_size(1) * out_size(2), 1);
    y = zeros(out_size(1) * out_size(2), 1);
    H_inv = inv(H);

    transformedImage = zeros(out_size);

    for j = 1:out_size(2)
        for i = 1:out_size(1)
            p = H_inv * [i+x_min; j+y_min; 1];
            a = round(p(1)/p(3));
            b = round(p(2)/p(3));
            if a>= 1 && b>= 1 && a<=s(1) && b<=s(2)
                transformedImage(i, j, :) = image(a,b,:);
            end
        end
    end

end
