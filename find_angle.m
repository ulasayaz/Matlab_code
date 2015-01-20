function [rad]=find_angle(x,y);

%finds angle between x and y in degrees (<90')
x=x/norm(x);y=y/norm(y);
rad=acos(x'*y);
% angle=180*rad/pi;
% if angle > 90
%     angle=180-angle;
% end
if rad > pi/2
    rad=pi-rad;
end