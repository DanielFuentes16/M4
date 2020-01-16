function [xhat] = triangulate( x1, x2, p1, p2, dim )

A = [x1(1)*p1(3,:) - p1(1,:);
    x1(2)*p1(3,:) - p1(2,:);
    x2(1)*p2(3,:) - p2(1,:);
    x2(2)*p2(3,:) - p2(2,:)];
[U, S, V] = svd(A);

xhat = V(:,4);
xhat = xhat./xhat(4);
end
