function [DTCmask] = DetectDTC(AbGPla, iminfo, pix, mupr)

%%% DTCmask: 
%%%        | 1: mask image on each z-plane 
DTCmask = cell(iminfo(4),2);
fint = cell(iminfo(4),1);
% totInt = zeros(iminfo(4),1);
mofm = zeros(iminfo(4),1);
% mupr = 1;

% make general threshold
parfor i = 1:iminfo(4)
    % calculate background for individual z planes and save in 'mofm'
    tt = reshape(AbGPla{i,1},[],1);
    tt(tt ==0) = [];
    mofm(i) = mean(tt) * 1;  
end

% remove NaN from 'mofm' to calculate 'mmofm'
mofm(isnan(mofm)) = 1;
mmofm = mean(mofm);

mofm = mofm * mupr;
mmofm = mmofm * mupr;

for i = 1:iminfo(4)
    t1 = AbGPla{i,1};
    t1(t1 < mofm(i)) = 0;

    t2 = medfilt2(t1, [3 3]);
    t2 = medfilt2(t2, [2 2]);

%     thrs = graythresh(t2)* 0.7;   %%% old setting for CHL
%     thrs = graythresh(t2)* 0.2;    %%% new setting for Sarah C
    thrs = graythresh(t2)* 10;    %%% for Sarah C Lng images
    
    if thrs > 1
        thrs = 1;
    end

    t3 = imbinarize(t2, thrs);
    t4 = bwareaopen(t3, round(pix*1.5));
    t5 = imclose(t4, true(pix*2));

    DTCmask{i,1} = t5;
    fint{i} = (AbGPla{i,1} - mofm(i)/mupr) / mofm(i)*mupr * mmofm;
%     DTCmask{i,2} = fint{i} .* uint8(t5);
    DTCmask{i,2} = fint{i} .* uint16(t5);
%     totInt(i) = sum(sum(uint8(t5) .* fint{i}));
end


%%% visaul DTC
% for i = 1:iminfo(4)
%     figure, imshow(DTCmask{i});
%     pause
%     close all
%     fprintf('\n%d-th zplane.', i);
% end


%%% reconstruct DTC in 3D.
% blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
%       | 6th: matching nuc #(row # in 'nuc')
% tfDTC = cat(3,DTCmask{:});
% temp = bwconncomp(tfDTC, 26);
% connDTC = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
% blDTC = struct2cell(connDTC)';

















