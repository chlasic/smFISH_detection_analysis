%%% gets signal intensity
function cirsum = GetIntensity(image, centers, radii) 

centers = round(centers);
radii = round(radii);
[xx, yy] = meshgrid(-radii:radii);
temp1 = sqrt(xx.^2 + yy.^2) <= radii;
temp2 = sqrt(xx.^2 + yy.^2) >= radii *5/10;
temp = temp1 .* temp2;
clear temp1 temp2

if centers(2) - radii < 1
    xa = 1;
    temp(:,1:radii-centers(2)+1) = [];
else xa = centers(2)-radii;
end
if centers(2) + radii > length(image(:,1))
    xb = length(image(:,1));
    temp(centers(2)+radii-length(image(:,1)):end,:) = [];
else xb = centers(2)+radii;
end
if centers(1) - radii < 1
    ya = 1;
    temp(:,1:radii-centers(1)+1) = [];
else ya = centers(1)-radii;
end
if centers(1) + radii > length(image(1,:))
    yb = length(image(1,:));
    temp(centers(1)+radii-length(image(1,:)):end,:) = [];
else yb = centers(1)+radii;
end

% if xa > xb
%     t1 = xb:xa;
% else
%     t1 = xa:xb;
% end
% 
% if ya > yb
%     t2 = yb:ya;
% else
%     t2 = ya:yb;
% end

% crop = image(t1, t2);
crop = image(xa:xb, ya:yb);
coor = temp == 1;
cirsum = sum(crop(coor));
% [ix, iy] = find(temp == 1);
% icoor = [ix iy] - radii - 1;
% icoor(:,1) = icoor(:,1) + centerss(1);
% icoor(:,2) = icoor(:,2) + centerss(2);
% 
% tsum = 0;
% for i = 1:length(icoor(:,1)
%     tsum = tsum + image(icoor(i,1), icoor(i,2));
% end
