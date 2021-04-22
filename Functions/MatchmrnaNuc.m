function [mrnas, mrnasNoATS, rnas_atsEx, nucm, nuc_atsEx] = MatchmrnaNuc(DETrna, thrNuc, nucs, iminfo, radius, vLim, pixr)
%%% This function matches rnas spots with nuclei (only mRNA spots)

%%% determine which nuc the RNA spots belong. 
% 'rnas' | 1:3 col = x,y,z coor | 4 col = associated nuc | 5th: vol ratio to mean vol per spot |
%        | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
%        | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
%        | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

mrnas = DETrna;
rnas_atsEx_all = mrnas;
mrnas(:,9:12) = 0;
nucm = nucs;
nuc_atsEx = nucm;
nucm(:,8) = 1:size(nucm,1);
nucm(:,9:11) = 0;
nuc_atsEx(:,8:9) = 0;
rnasNucList = cell(size(mrnas,1),1);
sen = strel('disk', pixr+1); 
meanIntensity_mRNA = mean(mrnas(:,6));  

% 'nuc' : | 1:3 xyz-coor | 4: radius | 5: ave. circularity |
%         | 6: std. circularity | 7: total DAPI sig. area (# pixel) | 8: nuc ID.
%         | 9: Voronoi mRNA count | 10: Voronoi count w/ limit (vLim, 3um) | 11: # mRNA in ROI w/ radius (2.5um)  |
idxForATS = zeros(9999,1);
locr = 1;
for i = 1:size(mrnas,1)
    cutset = nucm(nucm(:,3) >= mrnas(i,3) - radius/iminfo(6) & nucm(:,3) <= mrnas(i,3) + radius/iminfo(6),:);
    distemp = zeros(size(cutset,1),2);
    distemp(:,1) = sqrt( (cutset(:,1)-mrnas(i,1)).^2 + (cutset(:,2)-mrnas(i,2)).^2 + (cutset(:,3)-mrnas(i,3)).^2 );
    distemp(:,2) = cutset(:,8);
    distemp(:,3:6) = cutset(:,1:4);
    
    distemp = sortrows(distemp, 1);
    
    if ~isempty(distemp)
        %%% Test if detected exonal ATS spot is on DAPI (overlap).
        if ismember(round(mrnas(i,7)), [1 2 3 4]) == 1
            tb = cat(3,thrNuc{1:round(mrnas(i,7))+4,1});
        elseif ismember(round(mrnas(i,7)), [iminfo(4)-3 iminfo(4)-2 iminfo(4)-1 iminfo(4)]) == 1
            tb = cat(3,thrNuc{round(mrnas(i,7))-4:end,1});
        else
            tb = cat(3,thrNuc{round(mrnas(i,7))-4:round(mrnas(i,7))+4,1});
        end
        tb = max(tb,[],3);
        tx = imdilate(tb,sen);

        tstart = round(mrnas(i,2) - pixr*2-1):round(mrnas(i,2) + pixr*2+1);
        tend = round(mrnas(i,1) - pixr*2-1):round(mrnas(i,1) + pixr*2+1);

        tstart(tstart < 1) = [];
        tstart(tstart > iminfo(3)) = [];
        tend(tend < 1) = [];
        tend(tend > iminfo(2)) = [];

        txpart = tx(tstart, tend);
        
        %%% Check ATS & DAPI overlay
        dapiOL = 0;
        if sum(sum(txpart)) >= pixr*   2    % total # pixels
            dapiOL = 1;
        end
        
        %%% record exonal ATS seperately.
        % 'distemp': | 1: distance Nuc-RNA spot | 2: nuc ID | 
        %            | 3-5: xyz nuc coor | 6: nuc radius in pixels |
        % Select 3 shortest length from nuc center.
        if distemp(1,1) < distemp (1,6) * 0.95 && dapiOL == 1 
            if rnas_atsEx_all(i,6) > meanIntensity_mRNA  *  1.5  % cutoff for ATS (e.g. 1 for ATS with its intensity equal to 1 mRNA.
                rnas_atsEx_all(i,4) = distemp(1,2);
                rnas_atsEx_all(i,9) = rnas_atsEx_all(i,6) / meanIntensity_mRNA;
                
                %%% Delete ATS from mRNA list
                idxForATS(locr) = i;
                locr = locr + 1;
                nuc_atsEx(distemp(1,2),8) = nuc_atsEx(distemp(1,2),8) + rnas_atsEx_all(i,5);
                if nuc_atsEx(distemp(1,2),8) > 4
                    nuc_atsEx(distemp(1,2),8) = 4;
                end
                nuc_atsEx(distemp(1,2),9) = nuc_atsEx(distemp(1,2),9) + rnas_atsEx_all(i,6);
                nuc_atsEx(distemp(1,2),10) = nuc_atsEx(distemp(1,2),9) / meanIntensity_mRNA;
            end
        else
            %%% count mRNAs per cell (ID-ed by nucleus) with 3 diff ways.
            mrnas(i,4) = 1;
            % simple Voronoi count
            nucm(distemp(1,2),9) = nucm(distemp(1,2),9) + 1;
            mrnas(i,9) = 1;
            mrnas(i,12) = distemp(1,1); % unit: pixel

            % Voronoi w/ distance limit
            if distemp(1,1) <= vLim/iminfo(6)
                nucm(distemp(1,2),10) = nucm(distemp(1,2),10) + 1;
                mrnas(i,10) = 1;
            end

            % count # mRNA in ROI (rad: Radius; 2.5 um)
            valDist = distemp(distemp(:,1) < radius/iminfo(6),:);
            mrnas(i,11) = length(valDist(:,1));
            rnasNucList{i} = [rnasNucList{i} valDist(:,2)'];

            if ~isempty(valDist)
                nucm(valDist(:,2),11) = nucm(valDist(:,2),11) + 1;
            end
        end
    end
end
idxForATS(locr:end) = [];

rnas_atsEx = rnas_atsEx_all(idxForATS,:);
mrnasNoATS = mrnas;
mrnas(idxForATS,:) = [];

% keep calculating # mRNA in ROI.
for i = 1:size(mrnas,1)
    if ~isempty(rnasNucList{i})
        nucm(rnasNucList{i},11) = nucm(rnasNucList{i},11)-1+1/length(rnasNucList{i});
    end
end
        
  
