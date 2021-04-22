function [nucs, DTCnuc] = DetectNucBlobs(thrNuc, aftNuc, Nuc, iminfo, tnran)

%%% match detected sphere with bwconncomp-processed blobs to catch missed
% non-circular objects
for i=1:size(thrNuc,1)
    if isempty(thrNuc{i})
        thrNuc{i} = zeros(iminfo(3),iminfo(2));
    end
end

tfNuc = cat(3,thrNuc{:,1});
temp = bwconncomp(tfNuc, 26);
conNuc = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');

% blNuc: |Area | Centroid | BoundingBox | Image (stack) | z-projection | radius | circularity
%       the biggest radius is caluculated based on circle detection
blNuc = struct2cell(conNuc)';
[~,~,temp] = cellfun(@size, blNuc(:,4));
blNuc(temp > 20 | temp < 2, :) = [];

% ========= put z-proj of 3 slides (1st method/3) ======================
temp(temp > 20 | temp < 2, :) = [];
temp = [ceil(temp/6) round(temp/3) round(temp/2) round(temp/3*2) round(temp/6*5)];
for i = 1:size(temp,1)
    blNuc{i,5} = sum(blNuc{i,4}(:,:,temp(i,:)),3);
end   
blNuc(:,5) = cellfun(@(x) x > 0, blNuc(:,5), 'uni', 0);
%=======================================================================

% ====== % put mid slide (2nd method/3) ===========================
%     temp(temp > 7/iminfo(7) | temp < 1/iminfo(7), :) = [];
%     temp = round(temp/2);
%     for i = 1:length(temp)
%         blNuc{i,5} = blNuc{i,4}(:,:,temp(i));
%     end
%==================================================================

% =======% put z-projection (3rd method/3) =======================
%     blNuc(:,5) = cellfun(@(x) sum(x,3), blNuc(:,4), 'uni', 0); 
%     blNuc(:,5) = cellfun(@(x) x > 0, blNuc(:,5), 'uni', 0);
%=================================================================    

% circle detection on z-projected images (or mid slide)
temp = cell(length(blNuc),3);
tsiz = zeros(length(blNuc),1);
blNuc(:,8) = {0};
for i = 1:size(blNuc,1)
    [temp{i,1}, temp{i,2}, temp{i,3}] = imfindcircles(blNuc{i,5}, tnran, 'sensitivity', 0.98);
    if isempty(temp{i,1}) == 1
        temp{i,1} = 0;
        temp{i,2} = 0;
        temp{i,3} = 0;
    end

    if length(temp{i,2}) > 1
        tmp = find(temp{i,2} == max(temp{i,2}));
        temp{i,1} = temp{i,1}(tmp(1),:);
        temp{i,2} = temp{i,2}(tmp(1),:);
        temp{i,3} = temp{i,3}(tmp(1),:);
    end

    % mark blobs too small or too big
    tsiz(i) = mean(blNuc{i,3}(4:5));
    if tsiz(i) < 3 /iminfo(6)*0.65 ...  % 3um is mean size of nuc
            || tsiz(i) > 3 /iminfo(6)*1.5
        blNuc{i,8} = 99;
    end

%         if blNuc{i,3}(4) < 3 /iminfo(6)*0.5 || blNuc{i,3}(5) < 3 /iminfo(6)*0.5
%             blNuc{i,8} = 99;
%         end
end

blNuc(:,6:7) = temp(:,2:3);

% remove blobs not detected as circle (very loose detection)
temp = cell2mat(blNuc(:,6));
blNuc(temp == 0,:) = [];

% remove blobs too small or too big
temp = cellfun(@(x) x == 99, blNuc(:,8));
blNuc(temp,:) = [];

% visualize blob detected
%     for i=1:length(blNuc)
%         imshow(blNuc{i,5})
%         fprintf('\n%d', i);
%         pause
%     end


% 'blobNuc' : | x-coor | y | z | radius | ave. circularity |
%           | std. circularity | total DAPI sig. area (# pixel) | 
% of all nuclei
blobNuc = cell2mat(blNuc(:,2));
blobNuc(:,3) = blobNuc(:,3) * iminfo(7) / iminfo(6); % make voxel a cube

blobNuc(:,4) = cell2mat(blNuc(:,6));
temp = cell2mat(blNuc(:,3));
temp = temp(:,4:5);
temp = max(temp,[],2);
tmp = find(blobNuc(:,4) == 0);
blobNuc(tmp,4) = temp(tmp);

blobNuc(:,5) = cell2mat(blNuc(:,7));
blobNuc(:,6) = 0;
blobNuc(:,7) = cell2mat(blNuc(:,1));

temp = 1.5/iminfo(6); % 1.5 um is roughly mean radius for nuclei
tmp = blobNuc(:,4) < temp*0.85 | blobNuc(:,4) > temp*1.8;
blobNuc(tmp,:) = [];
blNuc(tmp,:) = [];

temp = mean(blobNuc(:,7)); %using mean # pixels
tmp = blobNuc(:,7) < temp*0.4;
blobNuc(tmp,:) = [];
blNuc(tmp,:) = [];

temp = mean(blobNuc(:,4)); % using mean radius
tmp = blobNuc(:,4) < temp*0.7 | blobNuc(:,4) > temp*1.5;
blobNuc(tmp,:) = [];
blNuc(tmp,:) = [];

blobNuc(:,8) = 999; % mark and find objects not overlapping with detected circles

% mark '0' when find overlapping 
for i=1:length(blobNuc(:,1))
    for j=1:length(Nuc(:,1))
        dCen = sqrt( (blobNuc(i,1)-Nuc(j,1))^2 + (blobNuc(i,2)-Nuc(j,2))^2 + (blobNuc(i,3)-Nuc(j,3))^2 );
        if blobNuc(i,4)/Nuc(j,4) > 0.5 && blobNuc(i,4)/Nuc(j,4) < 2
            if dCen < (blobNuc(i,4) + Nuc(j,4))*0.7
                blobNuc(i,8) = 0;
%                     fprintf('\na %d %d', i,j);
            end
        else
            temp = max(blobNuc(i,4), Nuc(j,4));
            if dCen < temp* 0.7
                blobNuc(i,8) = 0;
%                     fprintf('\nb %d %d', i,j);
            end
        end 

        subt = bwconncomp(blNuc{i,5}, 18);
        subc = regionprops(subt, 'Area');
        if length(subc(:,1)) > 1
            blobNuc(i,8) = 0;
        end
    end
end

%     imshow(blNuc{26,5})

% add nuclei not overlapping with circle detection.
coor = blobNuc(:,8) == 999;
temp = blNuc(coor,:);

iint = zeros(1,length(temp(:,1)));
parfor i=1:length(temp(:,1))
    iint(i) = 0;
    for j=1:temp{i,3}(6)
        cutimg = aftNuc{round(temp{i,3}(3)),1};
        cutimg = cutimg(round(temp{i,3}(2)):round(temp{i,3}(2))+temp{i,3}(5)-1, round(temp{i,3}(1)): round(temp{i,3}(1))+temp{i,3}(4)-1);

        iint(i) = iint(i) + sum(sum(temp{i,4}(:,:,j) .* double(cutimg)));
    end
end

adden = cell2mat(temp(:,2));

if ~isempty(adden)
    adden(:,3) = adden(:,3) * iminfo(7) / iminfo(6); % make voxel a cube
    adden(:,4) = cell2mat(temp(:,6));
    adden(:,5) = cell2mat(temp(:,7));
    adden(:,6) = 0;
    adden(:,7) = iint';

    % remove huge but little DNA-containing nuc
    temp = mean(Nuc(:,7)./Nuc(:,4)./1000);
    te2 = adden(:,7)./adden(:,4)./1000;
    adden(te2 < temp*0.3,:) = [];

    % final nuc information.
    nucs = [Nuc; adden];

else
    nucs = Nuc;
end
%%% 'nucs' labels : | col 1: x-coor | y | z | radius | col 5: ave. circularity |
%%%                 | col 6: std. circularity | col 7: total DAPI sig. area (# pixel) | 
    
%     cor = corr(nuc(:,4),nuc(:,7));
%     
%     scatter(nuc(:,4)*iminfo(6),nuc(:,7));
%     xlabel('Diameter of Nucleus (\mum)', 'fontsize',15);
%     ylabel('DAPI signal intensity (a.u.)', 'fontsize',15); 
%     text(1.4,650000, strcat('p = ',num2str(cor)));
%     hold on
%     scatter(adden(:,4),adden(:,7), 'r', '.');
    

% Remove double-detected nuclei (spheres overlapping in great portion)
if ~isempty(nucs)
    nucs(:,8) = 0;
    nCount = 1;
    for i=1:size(nucs,1)
        if nucs(i,8) == 0
            nucs(i,8) = nCount;
            nCount = nCount + 1;
            for j=i+1:size(nucs,1)
                if nucs(j,8) == 0
                    Cendist = sqrt( (nucs(i,1) - nucs(j,1))^2 + (nucs(i,2) - nucs(j,2))^2 ); 
                    Zdist = sqrt((nucs(i,3) - nucs(j,3))^2 );
%                     dist = sqrt( (nucs(i,1) - nucs(j,1))^2 + (nucs(i,2) - nucs(j,2))^2 + (nucs(i,3) - nucs(j,3))^2 ); 
                    if Cendist < 1.5/iminfo(6) && Zdist < 3/iminfo(6)
                        nucs(j,8) = nucs(i,8);
                    end
                end
            end
        end
    end
end
    
for i=1:max(nucs(:,8))
    if sum(nucs(:,8) == i) > 1
        maxR = max(nucs(nucs(:,8) == i,4));
        nucs(nucs(:,8) == i & nucs(:,4) ~= maxR,:) = [];
    end
end
    
nucs(:,8) = [];

%%% Record DTC location and remove it from the nuc list
%%% This section uses nuclei location and circularity to determine the DTC.
%%% Comment out below if DTC detection does not work
nucs_sorted = sortrows(nucs, [1 5 6]);
if nucs_sorted(2,5)/nucs_sorted(1,5) < 0.5
    DTCnuc = nucs_sorted(2,:);
    dcs = 1;
elseif nucs_sorted(2,5)/nucs_sorted(1,5) > 2
    DTCnuc = nucs_sorted(1,:);
    dcs = 1;
else
    DTCnuc = [];
    dcs = 0;
end

if dcs == 1
    idx = ismember(nucs, DTCnuc, 'rows');
    nucs(idx,:) = [];
end



       
