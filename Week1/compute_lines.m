function [l1, l2, l3, l4] = compute_lines(p1, p2, p3, p4, p5, p6, p7, p8)
%COMPUTE_LINES Summary of this function goes here
%   Detailed explanation goes here

    point1 = [p1(1) p1(2) 1];
    point2 = [p2(1) p2(2) 1];
    point3 = [p3(1) p3(2) 1];
    point4 = [p4(1) p4(2) 1];
    point5 = [p5(1) p5(2) 1];
    point6 = [p6(1) p6(2) 1];
    point7 = [p7(1) p7(2) 1];
    point8 = [p8(1) p8(2) 1];

    l1 = cross(point1, point2);
    l2 = cross(point3, point4);
    l3 = cross(point5, point6);
    l4 = cross(point7, point8);

end

