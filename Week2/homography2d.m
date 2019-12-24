function [H] = homography2d(x1, x2)
% to compute the matrix H, we should construct matrix A
A=[];
tempA=[];
% for every point, we have only two linearly independent equations
% I combine them together and get a 12x9 matrix A
[p1_norm, p1_trans]= normalise2dpts(x1);
[p2_norm, p2_trans]= normalise2dpts(x2);
for i=1:4
    %tempA=[zeros(1,3) (-1)*x2(:,i)' x1(2,i)*x2(:,i)';
    %       x2(:,i)' zeros(1,3) -x1(1,i)*x2(:,i)'];
    tempA=[zeros(1,3) (-1)*p2_norm(:,i)' p1_norm(2,i)*p2_norm(:,i)';
        p2_norm(:,i)' zeros(1,3) -p1_norm(1,i)*p2_norm(:,i)'];
    A=[A ;tempA];
end

% conduct SVD on matrix A
[U,S,V]=svd(A);
% the right-most-column of V is the h with least squares
h=V(:,end);
H=reshape(h,[3, 3])';
% denormalise H, according to the transformation matrixs
H = inv(p1_trans)*H*(p2_trans);

end