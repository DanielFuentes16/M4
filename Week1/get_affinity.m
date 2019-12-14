function [H] = get_affinity(l1, l2, l3, l4)
%AFFINITY Summary of this function goes here
%   Detailed explanation goes here

vanishing_point1 = cross(l1,l2);
vanishing_point1 = vanishing_point1/vanishing_point1(3);
vanishing_point2 = cross(l3,l4);
vanishing_point2 = vanishing_point2/ vanishing_point2(3);
vanishing_line = cross(vanishing_point1,vanishing_point2);

H = [1 0 0;
     0 1 0; 
     vanishing_line(1)/vanishing_line(3) vanishing_line(2)/vanishing_line(3) 1];
end