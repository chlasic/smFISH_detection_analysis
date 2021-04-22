function [blobrna, blrna, mom] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, numOfImg, g, thresForIntron)
% This function detects intronic RNA spots (transcription sites). The
% output is RNA coordinates, counts and intensity

sizel = iminfo(4); 
t1 = cell(sizel,1); % cropped RNA images only in the nuclei.
mofm = zeros(sizel,1);
thrrna = cell( sizel, 1);
Pmofm = zeros(sizel,1);
sen = strel('disk', pixr*8); 
se = strel('disk', pixr*1); 
se2 = strel('disk', pixr*2); 
mupr = thresForIntron + 0; % multiplier to define the background level (becomes higher when bg is high)
if thresForIntron <= 0
    mupr = 0.001;
end

% make general threshold
parfor i = 1:sizel
    t1{i} = rPla{i,1};
    
    % black out RNA images not overlapping with nuclei
    if ismember(i, [1 2 3]) == 1
        tb = cat(3,thrNuc{1:i+3,2})*1.5;
    elseif ismember(i, [iminfo(4)-2 iminfo(4)-1 iminfo(4)]) == 1
        tb = cat(3,thrNuc{i-3:end,2})*1.5;
    else
        tb = cat(3,thrNuc{i-3:i+3,2});
    end
    tb = max(tb,[],3);
    tx = imdilate(tb,sen);
    
    if isa(t1{i},'uint16')
        modr = uint16(tx);
    elseif isa(t1{i},'uint8')
        modr = uint8(tx);
    end
    p1{i} = t1{i} .* modr;
    
    % make a binary image where background signal is.
    thrs = graythresh(p1{i}(p1{i}>0));
    if thrs == 0
        thrs = 0.001;
    end
    mask1 = imbinarize(p1{i}, thrs );
    
    % calculate background level only inside the gonad using gonadal mask.
    % calculate background for individual z planes and save in 'mofm'
    tAll = p1{i}(mask1 > 0);
    mofm(i) = mean(tAll) * 1.02;    
end

% remove NaN from 'mofm' to calculate 'mmofm'
tempM = mean(mofm(~isnan(mofm)));
mofm(isnan(mofm)) = tempM;
mofm = mofm * mupr;
mom = mean(mofm);


%%% threshold images and detect RNA spots
parfor i = 1:sizel
%%%----------------- peak detection method (1/3) ------------------
    [A, B] = FastPeakFind(t1{i}, mofm(i));
    A = [A(1:2:end) A(2:2:end)];

    if ~isempty(A)
        %%% Overall signal should be 25% brighter than the overall background
        %%% This also looks at individual detected peaks and compare to its surrounding.
        %%% That ratio should be larger than 1.1
        temp = mean(t1{i}(B == 1));
        if isnan(temp)
            Pmofm(i) = 0;
        else
            Pmofm(i) = temp;
        end

        A(:,3) = 0;
        if Pmofm(i)/mofm(i) < 1.05 * mupr   % cutoff for overall dots-bg ratio
            B(A(:,2), A(:,1)) = 0;
            A(:,3) = 999;
        else
            ratioSB = zeros(size(A,1),2);
            for j = 1:size(A,1)
                rng1 = pixr*2;
                rng2 = pixr*10;
                
                Xrange1 = A(j,1)-rng1:A(j,1)+rng1 ;
                Yrange1 = A(j,2)-rng1:A(j,2)+rng1;
                Xrange1 = Xrange1(Xrange1 > 0 & Xrange1 <= iminfo(2));
                Yrange1 = Yrange1(Yrange1 > 0 & Yrange1 <= iminfo(3));
                sigMean1 = t1{i}( Yrange1 , Xrange1 );
                sigMean2 = mean(mean(sigMean1));
                
                Xrange2 = A(j,1)-rng2:A(j,1)+rng2;
                Yrange2 = A(j,2)-rng2:A(j,2)+rng2;
                Xrange2 = Xrange2(Xrange2 > 0 & Xrange2 <= iminfo(2));
                Yrange2 = Yrange2(Yrange2 > 0 & Yrange2 <= iminfo(3));
                bgMean1 = t1{i}( Yrange2 , Xrange2 );
                bgMean1(rng2-rng1*2+1:rng1+rng2+1, rng2-rng1*2+1:rng1+rng2+1) = 0;
                bgMean2 = bgMean1(bgMean1 > 0);
                sortN2 = sort(bgMean2(:));
                bgMean3 = mean(sortN2(10:end));
                ratioSB(j,1) = sigMean2/bgMean3;
%                 fprintf('\n%.2f for %d', ratioSB(j), j)
                
                if ratioSB(j,1) < 1.15 * mupr  % the cutoff signal-local bg ratio of each ATS (3x3 ATS vs. 7x7 bg region).
                    B(A(j,2), A(j,1)) = 0;
                    A(j,3) = 999;
                    ratioSB(j,2) = 999;
                end
            end
        end
        A(A(:,3) == 999,:) = [];
        ratioSB(ratioSB(:,2) == 999,:) = [];
        
        %%% Quick visual -----------------------
%         figure, imshow(t1{i}*100)
%         hold on
%         plot(A(:,1), A(:,2), 'r+');
%         for j = 1:size(A,1)
%             text(A(j,1), A(j,2), num2str(ratioSB(j)), 'color', 'c');
%         end
        %%%-------------------------------------
      
        if ~isempty(A)
            A(:,3) = [];

            % remove peaks too close
            A(:,3) = 0;
            for j = 1:length(A(:,1))
                for k = j+1:length(A(:,1))
                    dist = sqrt((A(k,1)-A(j,1))^2 + (A(k,2)-A(j,2))^2);
                    if dist < pixr*2
                        A(j,3) = 999;
                        B(A(j,2),A(j,1)) = 0;
                    end
                end
            end
            A(A(:,3) == 999,:) = [];
            A(:,3) = [];

%%%---------- Image segmentation (bwconncomp) method (2/3) ---------
            t2 = t1{i} - mofm(i) * 0.8;
            t2(t2 < 0 ) = 0;

            loopn = 2;
            for j = 1:loopn
                t2 = medfilt2(t2, [pixr pixr]);
                tt = reshape(t2,[],1);
                tt = sort(tt);
                ty = tt(tt ~= 0);
                tmofm = mean(mean(ty(20:end)));
                t2(t2 < tmofm) = 0;
            end
% 
%             tpix = t2(end,:);
%             t2 = t2(1:end-1,:);
%             t2 = [tpix; t2];
%             tpix = t2(:,end);
%             t2 = t2(:,1:end-1);
%             t2 = [tpix t2];

            thrs = graythresh(t2);
            if thrs == 0
                thrs = 0.001;
            end

            t3 = imbinarize(t2, thrs );
            t5 = bwareaopen(t3, ceil(pixr*2));
%             t5 = t4 & ~bwareaopen(t4, ceil(pixr*60));
            t6 = imdilate(t5, se2);

%%%------ Cross validate detected spots (bwconncomp + peakfinder) ----
%%% This removes bwconncomp-detected dots that are not found by peakfinder.
            A(:,3) = 0;
            for j=1:size(A,1)
                if t6(A(j,2),A(j,1)) == 0
                    A(j,3) = 999;
                    B(A(j,2),A(j,1)) = 0;
                end
            end
            A(A(:,3) == 999,:) = [];
            
            if ~isempty(A)
                A(:,3) = [];
                
                
%=========== visual detected RNAs =======
%             figure,imshow(t1{i}*100);
%             hold on
%             plot(A(:,1), A(:,2), 'r+');
% ---------------------------------------
        
        

            %%%-------- signal-to-local background ratio method (3/3) ----------
            % get intensity of spots and background to get the s:n ratio.
                sizC = pixr*2-1; % radius of the center region to get intensity
                cal = zeros(size(A,1), 4);
                if ~isempty(A) && sum(sum(A)) ~= 0
                    for winSize = 2:4
                        for j=1:size(A,1)
                            % 'A' col: | x-coor | y-coor | center_intensity | backgd2 int | bkgd3 int | bkgd4 int |.
                            sizB = pixr*2*winSize; 
                            xb = 1:sizB*2+1;

                            [x, y] = meshgrid(xb,xb);
                            Rb = sqrt((x-sizB-1).^2 + (y-sizB-1).^2); % center circle

                            Rc = Rb;
                            Rc(Rc>sizC) = 0;
                            Rc(sizB+1,sizB+1) = 1;
                            Rb(Rb>sizB) = 0;
                            Rb(Rb<sizC+2) = 0;

                            LB = A(j,2) - sizB; % left boundary of the ROI for scanning
                            if LB < 1
                                Rb(:,1:1-LB)=[];
                                Rc(:,1:1-LB)=[];
                                LB = 1;
                            end

                            RB = A(j,2) + sizB;
                            if RB > iminfo(3)
                                Rb(:,end-RB+iminfo(3)+1:end)=[];
                                Rc(:,end-RB+iminfo(3)+1:end)=[];
                                RB = iminfo(3);
                            end

                            UB = A(j,1) - sizB;
                            if UB < 1
                                Rb(1:1-UB,:)=[];
                                Rc(1:1-UB,:)=[];
                                UB = 1;
                            end

                            DB = A(j,1) + sizB;
                            if DB > iminfo(2)
                                Rb(end-DB+iminfo(2)+1:end,:)=[];
                                Rc(end-DB+iminfo(2)+1:end,:)=[];
                                DB = iminfo(2);
                            end

                            cutcen = t1{i}(LB:RB,UB:DB);
                            centr = cutcen(Rc > 0);

                            cutbck = t1{i}(LB:RB,UB:DB);
                            surrd = cutbck(Rb > 0);

                            % subtract detected spots that are dim (not significant spots).
                            [~,cal(j,1)] = kstest2(centr, surrd);
                            cal(j,2) = mean(centr)/mean(surrd); %ratio to surrounding bg.
                            cal(j,3) = mean(centr)/mofm(i); % ratio to mean background (fixed)
                            cal(j,4) =  j;
                            caller = 1;

                            if winSize == 2
                                if cal(j,1) > 5*10^-1 / mupr || cal(j,2) < 1.05 * mupr || cal(j,3) < 1.05 * mupr
                                    B(A(j,2),A(j,1)) = 0;
                                    caller = 0;
                                end
                            elseif winSize == 3
                                if cal(j,1) > 1*10^-1 / mupr || cal(j,2) < 1.1 * mupr
                                    B(A(j,2),A(j,1)) = 0;
                                    caller = 0;
                                end
                            elseif winSize == 4
                                if cal(j,1) > 5*10^-2 / mupr || cal(j,2) < 1.15 * mupr
                                    B(A(j,2),A(j,1)) = 0;
                                    caller = 0;
                                end
                            end

            % =========== show results ========================
%                             fprintf('\n%dround %dth-> p: [[%d]] %d, ratio: %d, ratio(pos back): %d', winSize, j, caller, cal(j,1), cal(j,2), cal(j,3));
%                             if j == size(A,1)
%                                 fprintf('\n   < %dXradius applied >',winSize);
%                                 fprintf('\n');
%                             end
            %---------- save results for review -----------
                        end
                    end
                end
                thrrna{i} = imdilate(B, se);
            end
        end
    end
    fprintf('\n%d(th) Ch: %d(th)/total %d images,... %d(th)/ %d z-planes.', g, f, numOfImg, i, sizel);
end
fprintf('\n');



%%% 3D reconstitution using bwconncomp.
% blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
%        | 6th: matching nuc #(row # in 'nuc')
for i=1:size(thrrna,1)
    if isempty(thrrna{i})
        thrrna{i} = zeros(iminfo(3),iminfo(2));
    end
end
        
tfrna = cat(3,thrrna{:});
temp = bwconncomp(tfrna, 26);
conrna = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
blrna = struct2cell(conrna)';

% calculates intensity of detected blobs from original images.
iint = zeros(1,length(blrna(:,1)));
for i=1:length(blrna(:,1))
    iint(i) = 0;
    rng = round( blrna{i,3}(3) - (blrna{i,3}(6)-1)/2:blrna{i,3}(3) + (blrna{i,3}(6)-1)/2  );
    rng = rng(rng > 0 & rng <= iminfo(4));
    for j=1:length(rng)
        wholeImg = t1{rng(j)};
        rangeY = round(blrna{i,3}(2)):round(blrna{i,3}(2))+blrna{i,3}(5)-1;
        rangeX = round(blrna{i,3}(1)): round(blrna{i,3}(1))+blrna{i,3}(4)-1; 
        cutImg = wholeImg( rangeY  , rangeX  );
%         plot( rangeX(1) , rangeY(1), 'r+');
        iint(i) = iint(i) + sum(sum(blrna{i,4}(:,:,j) .* double(cutImg)));
    end
end

blrna(:,5) = num2cell(iint');
% fprintf('\nsize: %d %d\n', size(blrna,1), size(blrna,2));    

if  isempty(blrna) == 0
    if isempty(blrna{1,1}) == 0
        % 'blobrna': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
        %            | col 6: total intensity | x centroid | y centroid | z centroid | col 10: # ATS in the blob.
        blobrna = zeros(length(blrna(:,1)),1);
        blobrna(:,1) = 1:length(blrna(:,1));
        temp = cell2mat(blrna(:,3));
        blobrna(:,2:4) = temp(:,4:6);
        temp = cell2mat(blrna(:,2));
        blobrna(:,7:9) = temp(:,1:3);

        blobrna(:,6) = cell2mat(blrna(:,5));
        blobrna(:,5) = cell2mat(blrna(:,1));
  
        %%% remove false-positively detected blobs
        if isempty(blobrna) == 0
            % get mean of obj. size (pixels) and intensity (normalized: ind. int. - mean)
            mint = mean(blobrna(:,6));
%             mpix = mean(blobrna(:,5));

            % Exclude RNA spots seen only on 1 z-plane and dim.
            temp = blobrna(:,4) < 2 & blobrna(:,6) < mint * 0.2 ;
            blobrna(temp,:) = [];

            if isempty(blobrna) == 0
                %%% Exclude spots smaller than 2x2 pixels
                temp = blobrna(:,2) < pixr*2 & blobrna(:,3) < pixr*2;
                blobrna(temp ,:) = [];
            end
        end
    else
        blobrna = [];
    end
else
    blobrna = [];
end



%%% Quick visual
% figure, imshow(rZproj*50)
% hold on
% for i = 1:size(blobrna,1)
%     text(blobrna(i,7), blobrna(i,8), num2str(i), 'color', 'c');
% %     plot(blobrna(:,7), blobrna(:,8), 'r+');
% end







