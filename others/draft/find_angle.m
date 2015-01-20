function [rad]=find_angle(x,y);

%finds angle between x and y in radian ( < pi/2)
a=norm(x);b=norm(y);c=norm(x-y);
T=((a^2+b^2-c^2)/(2*a*b));
rad=acos(T);
% angle=180*rad/pi;
% if angle > 90
%     angle=180-angle;
% end
if rad > pi/2
    rad=pi-rad;
end