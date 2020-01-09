function [F] = fundamental_matrix(points1 ,points2)
%FUNDAMENTAL_MATRIX Calculates fundamental matrix using the 8-point algorithm
%  
    N = size(points1, 2);
    W = ones(N, 9);
    for i=1:N
        u = points1(1,i) / points1(3,i);
        v = points1(2,i) / points1(3,i);
        u_prime = points2(1,i) / points2(3,i);
        v_prime = points2(2,i) / points2(3,i);
        W(i, :) = [u*u_prime v*u_prime u_prime u*v_prime v*v_prime v_prime u v 1];
    end
    [U,D,V] = svd(W);
    f = V(:,9);
    Fr3 = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];
    [U,D,V] = svd(Fr3);
    D(3,3) = 0;
    F = U * D * V';

end

