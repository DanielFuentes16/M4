function [transformedImage] = apply_H(image, H)
% Apply a homography to an image
%

image = double(image);  
s = size(image);

topleft = H * [0; 0; 1];        
topright = H * [s(2); 0; 1];     
bottomleft = H * [0; s(1); 1];     
bottomright = H * [s(2); s(1); 1];  

topleft = topleft / topleft(3);
topright = topright / topright(3);
bottomleft = bottomleft / bottomleft(3);
bottomright = bottomright / bottomright(3);

x1 = max([topleft(1), topright(1), bottomleft(1), bottomright(1)]);
y1 = max([topleft(2), topright(2), bottomleft(2), bottomright(2)]);
out_size = round([y1, x1, 3])

x = zeros(out_size(1) * out_size(2), 1);
y = zeros(out_size(1) * out_size(2), 1);
H_inv = inv(H);

for j = 1:out_size(2)
    for i = 1:out_size(1)
        p = H_inv * [j; i; 1];
        idx = (out_size(1) * (j - 1)) + i;
        x(idx) = p(1) / p(3);
        y(idx) = p(2) / p(3);
    end
end

ch1 = interp2(image(:,:,1), x, y);
ch2 = interp2(image(:,:,2), x, y);
ch3 = interp2(image(:,:,3), x, y);

transformedImage = zeros(out_size);
transformedImage(:,:,1) = reshape(ch1, out_size(1), out_size(2), 1);
transformedImage(:,:,2) = reshape(ch2, out_size(1), out_size(2), 1);
transformedImage(:,:,3) = reshape(ch3, out_size(1), out_size(2), 1);

end

