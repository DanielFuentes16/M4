function angle = compute_angle(line1 ,line2)
%COMPUTE_ANGLE Summary of this function goes here
%   Detailed explanation goes here

angle = real(180 - acos((line1(1)*line2(1)+line1(2)*line2(2))/(sqrt(line1(1)^2+line2(1)^2)*sqrt(line1(2)^2+line2(2)^2)))*180/pi);

%m2 = line2(2)/line2(1);
%m1 = line1(2)/line1(1);
%angle = atan(abs((m2-m1)/(1+m1*m2)));

end