function [cyl_axis, t6] = DefOutline(PhaPla, pix, iminfo)
% This function finds cell outlines based on input phase images and outputs
% coordinates and radii (to the outline) of cell center line. Input
% variables are following. PhaPla: phase image stack (taken out from
% confocal images), pix: smallest pixel size control (automatically 
% determined), iminfo: image info imported from metafile (automatically
% done).

% se = strel('disk',2); 

% genenrate z-projected image (all planes)
% lon = 1:length(PhaPla(:,1));

% use 1/3 middle z-slices to make z-projected phase image
lon = round(length(PhaPla(:,1))/3):round(length(PhaPla(:,1))/3*2);
temp = PhaPla{lon,1};
zproj = uint8(sum(cat(3, temp),3)/length(PhaPla));


% apply image filters for crisper images
t1 = imsharpen(zproj, 'radius', pix*4, 'amount', pix*4); % sharpen
t2 = imcomplement(t1); % invert image
t3 = medfilt2(t2, [pix*2 pix*2]);
t4 = t3 == 255;
% t4 = imerode(t4,se);

t4(:,iminfo(2)-2:iminfo(2)) = cat(2, t4(:,iminfo(2)-3), t4(:,iminfo(2)-3), t4(:,iminfo(2)-3));

%     t5 = imfill(t4, [iminfo(2)/2 + 20 iminfo(3)/2],8);
%     t5 = imfill(t5, [iminfo(2)/2 - 20 iminfo(3)/2],8);
t6 = bwareaopen(t4, pix*500, 8);

outcoor = bwboundaries(t6);
outsiz = cell2mat(cellfun(@length, outcoor, 'uni', 0));
outsiz = outsiz(:,1);
% outlc contains x,y coor of germline outline.
outlc = cell2mat(outcoor(outsiz == max(outsiz)));
outlc(:,[1,2]) = outlc(:,[2,1]);

% reconsititute as 3D cylinder
outlc(:,3) = round(iminfo(4)/2);

%%% gets cylinderical axis coordinates and r to draw the cell on the plot.
%%% "cline" is a sorted Cell outline coordinates.
%%% (1) x-coor | (2) y-coor | (3) z-coor | (4) Nb mask points on that x |
%%%
%%% "cyl_axis" contains: (unit: pixels)
%%% (1) x-coor | (2) y-coor | (3) z-coor | (4) radius of the cell on that point.| 
cline = sortrows(outlc, 1);
cline(:,4) = 1;
cline = sortrows(cline, [1 2]);

for i=2:length(cline(:,1))
    if cline(i,1) == cline(i-1,1)
        cline(i,4) = cline(i,4) + cline(i-1,4);
    end
end
cline = sortrows(cline, [1 -4]);

lent = cline(end,1)-cline(1,1)+1; 
jr = 1;
cyl_axis= zeros(lent,4);
for i=1: lent
    cyl_axis(i,1:4) = [cline(jr,1) (cline(jr,2)+cline(jr+cline(jr,4)-1,2))/2 cline(jr,3) cline(jr,2)-cline(jr+cline(jr,4)-1,2)/2];
    jr = jr + cline(jr,4);
end