function [rnas, nucm] = MatchrnaNuc(thrNuc, blobrna, nucs, iminfo, pixr, dists)
% This function matches rnas spots with nuclei (only internuclear spots)

%%% determine which nuc the RNA spots belong. 
% 'rnas' | 1:3 col = x,y,z coor | 4 col = nuc #. | 5th: 1 or 2 spots |
%        | 6th: total sig intensity (trx activity) | 7: z-coor in plane #.
rnas = blobrna(:,7:9);
rnas(:,3) = rnas(:,3) * iminfo(7) / iminfo(6);
rnas(:,4:5) = 0;
rnas(:,5) = floor(blobrna(:,5) / mean(blobrna(:,5)) *   3/4   );
rnas(rnas(:,5) == 0,5) = 1;
rnas(:,6) = blobrna(:,6);
rnas(:,7) = blobrna(:,9);
% hotnuc = cell(length(nuc(:,1)),1);
nucm = nucs;
nucm(:,8:9) = 0;
sen = strel('disk', pixr+1); 

for i=1:size(rnas,1)
    distemp = zeros(size(nucm,1),2);
    for j=1:size(nucm,1)
        distemp(j,1) = sqrt( (nucm(j,1)-rnas(i,1))^2 + (nucm(j,2)-rnas(i,2))^2 + (nucm(j,3)-rnas(i,3))^2 );
        distemp(j,2) = j;
    end

    % 'distemp': | distNR | nuc ID | nuc radius |
    % select 3 shortest length from nuc center.
    distemp = sortrows(distemp, 1);
    if size(distemp,1) >=3
        distemp = distemp(1:3,:);
    end
    distemp(:,3) = nucm(distemp(:,2),4);

    distemp(distemp(:,1) > distemp(:,3) *    dists   ,:) = []; 
    
    
    % Test if detected ATS spot is on DAPI (overlap).
    if ismember(round(rnas(i,7)), [1 2 3 4]) == 1
        tb = cat(3,thrNuc{1:round(rnas(i,7))+4,1});
    elseif ismember(round(rnas(i,7)), [iminfo(4)-3 iminfo(4)-2 iminfo(4)-1 iminfo(4)]) == 1
        tb = cat(3,thrNuc{round(rnas(i,7))-4:end,1});
    else
        tb = cat(3,thrNuc{round(rnas(i,7))-4:round(rnas(i,7))+4,1});
    end
    tb = max(tb,[],3);
    tx = imdilate(tb,sen);
    
    tstart = round(rnas(i,2) - pixr*2-1):round(rnas(i,2) + pixr*2+1);
    tend = round(rnas(i,1) - pixr*2-1):round(rnas(i,1) + pixr*2+1);
    
    tstart(tstart < 1) = [];
    tstart(tstart > iminfo(3)) = [];
    tend(tend < 1) = [];
    tend(tend > iminfo(2)) = [];
    
    txpart = tx(tstart, tend);

    if isempty(distemp) == 1
        rnas(i,4) = 0;
    elseif sum(sum(txpart)) < pixr*  2 
        rnas(i,4) = 0;
    else
        rnas(i,4) = distemp(1,2);
        nucm(distemp(1,2),8) = nucm(distemp(1,2),8) + rnas(i,5);
        if nucm(distemp(1,2),8) > 4
            nucm(distemp(1,2),8) = 4;
        end
        nucm(distemp(1,2),9) = nucm(distemp(1,2),9) + rnas(i,6);
    end
end     


% remove rnas that are not in a nucleus
rnas(rnas(:,4) == 0,:) = [];
