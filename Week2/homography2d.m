function [H] = homography2d(x1, x2)
    % to compute the matrix H, we should construct matrix A
    x1 = normalise(x1);
    x2 = normalise(x2);

    % Compute a similarity transformation T1 & T2
    mu = mean(x1, 2);
    x1t = x1(1,:) - mu(1);
    x2t = x1(2,:) - mu(2);
    d = sqrt(x1t.^2 + x2t.^2);
    sigma = sqrt(2) / mean(d);
    T1 = [sigma, 0, -mu(1)*sigma; 0, sigma, -mu(2)*sigma; 0, 0, 1];
    
    mu = mean(x2, 2);
    x1t = x2(1,:) - mu(1);
    x2t = x2(2,:) - mu(2);
    d = sqrt(x1t.^2 + x2t.^2);
    sigma = sqrt(2) / mean(d);
    T2 = [sigma, 0, -mu(1)*sigma; 0, sigma, -mu(2)*sigma; 0, 0, 1];
    
    % The centroid of the transformed points T*x is the coordinate origin   
    xn1 = T1 * x1;
    xn2 = T2 * x2;

    A = [];
    for i = 1:size(xn1, 2)
        a = -xn2(3,i) * xn1(:,i);
        b = xn2(1,i) * xn1(:,i);
        c = xn2(2,i) * xn1(:,i);
        A = [A; a', 0, 0, 0, b'; 0, 0, 0, a', c'];
    end

    % conduct SVD on matrix A
    [~, ~, V] = svd(A);
    
    % the right-most-column of V is the h with least squares
    H_hat = reshape(V(:,end), 3, 3)';
    
    % denormalise H, according to the transformation matrixs
    H = inv(T2) * H_hat * T1;
end

function xn = normalise(x)
    xn = x ./ repmat(x(end,:), size(x,1), 1);
end