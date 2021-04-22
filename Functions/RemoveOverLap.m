function [centersNew,radiiNew,cirNew, intenNew]=RemoveOverLap(centers,radii, cir, inten, tol,option)
% This function deals with overlaping circles by:
% option 1: removes one circle of the two (if it does not matter which one).
% option 2: removes the smaller circle of the two
% option 3: kills both all overlaping circles
% option 4: kills one containing more black (or white) space. 
% centers - (x,y) circles centers.
% radii - the circles radius
% tol - tolerance for an overlap, im number of pixels.
% option - 1,2 or 3, 4, see above.
% Uses the function snip() from the file exchange.
l=length(centers);
if l >= 4
for i= 1: l
    s=i+1;
    for j=s:l
        d_ij=sqrt((centers(i,1)-centers(j,1)).^2+(centers(i,2)-centers(j,2)).^2);
        k=radii(i)+radii(j)-tol;
        if d_ij < k && radii(j)>0
            %option 1
            if option == 1
             centers(i,1)=0;
             centers(i,2)=0;
             radii(i)=0;
            end
            %option 2 
            if option == 2
                 if radii(i)>radii(j)
                     centers(j,1)=0;
                     centers(j,2)=0;
                     radii(j)=0;
                 else
                    centers(i,1)=0;
                    centers(i,2)=0;
                    radii(i)=0;
                 end
            end
            %option 3
            if option ==3
                     centers(j,1)=0;
                     centers(j,2)=0;
                     radii(j)=0;
                     centers(i,1)=0;
                     centers(i,2)=0;
                     radii(i)=0;
            end
            if option ==4
%                 inten(i,1) = GetIntensity(image, centers(i,:), radii(i));
%                 inten(j,1) = GetIntensity(image, centers(j,:), radii(j));
                if inten(i,1)/inten(j,1) > 0.83 && inten(i,1)/inten(j,1) < 1.17
                    if radii(i) >= radii(j)
                        centers(j,1) = 0;
                        centers(j,2) = 0;
                        radii(j) = 0;
                        cir(j) = 0;
                        inten(j) = 0;
                    else
                        centers(i,1)=0;
                        centers(i,2)=0;
                        radii(i)=0;
                        cir(i) = 0;
                        inten(i) = 0;
                    end
                
                elseif inten(i,1) > inten(j,1)
                    centers(j,1) = 0;
                    centers(j,2) = 0;
                    radii(j) = 0;
                    cir(j) = 0;
                    inten(j) = 0;
                else
%                 elseif inten(i,1) < inten(j,1) 
                     centers(i,1)=0;
                     centers(i,2)=0;
                     radii(i)=0;
                     cir(i) = 0;
                     inten(i) = 0;
                end
                
                % removes circle has low mean sig intensity.
%                 if inten(i) < 80
%                     centers(i,1)=0;
%                     centers(i,2)=0;
%                     radii(i)=0;
%                     cir(i) = 0;
%                     inten(i) = 0;
%                 end
%                 if inten(j) < 80
%                     centers(j,1) = 0;
%                     centers(j,2) = 0;
%                     radii(j) = 0;
%                     cir(j) = 0;
%                     inten(j) = 0;
%                 end
            end
%                 [xx, yy] = meshgrid(-radii(i):radii(i));
%                 temp = sqrt(xx.^2 + yy.^2) <= radii(i);
%                 [ix, iy] = find(temp == 1);
%                 icoor = [ix iy] - radii(i) - 1;
%                 icoor(:,1) = icoor(:,1) + centers(i,1);
%                 icoor(:,2) = icoor(:,2) + centers(i,2);
%                 inten(i,1) = sum(image(icoor(:,1), icoor(:,2)));
                
        end
   end
end
end
%create new circles vectors using snip()
centersNew=snip(centers,'0');
radiiNew=snip(radii,'0');
cirNew=snip(cir,'0');
intenNew=snip(inten,'0');
