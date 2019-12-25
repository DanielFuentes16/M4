function [H] = homography2d(x1, x2)
    % to compute the matrix H, we should construct matrix A
    A=[];

    % for every point, we have only two linearly independent equations
    % I combine them together and get a 12x9 matrix A
    [x1_norm, x1_t]= normalise2dpts(x1);
    [x2_norm, x2_t]= normalise2dpts(x2);

    for i=1:4
        tempA=[zeros(1,3) (-1)*x2_norm(:,i)' x1_norm(2,i)*x2_norm(:,i)';
               x2_norm(:,i)' zeros(1,3) -x1_norm(1,i)*x2_norm(:,i)'];
        A=[A ;tempA];
    end

    % conduct SVD on matrix A
    [U,S,V]=svd(A);

    % the right-most-column of V is the h with least squares
    h=V(:,end);
    H=[h(1) h(2) h(3);
       h(4) h(5) h(6);
       h(7) h(8) h(9);];

    % denormalise H, according to the transformation matrixs
    H = x1_t\H*(x2_t);
end