%%% This code detects smFISH dots and records their information in 3D, 
%%% including locations (xyz coordinates), intensities, information of nuclei. 
%%% This code is optimized for C. elegans germline but should work for any
%%% tissues by adjusting input parameters such as nuclear size range (variable 'nrange').
%%% The code reads tif (Tiff) images and lif (Leica) images. 
%%% This code is developed by ChangHwan Lee (SUNY Albany). Contact him at
%%% chlee@albany.edu for any questions.

%% Designate folder locations for image input and result storage.
clc; clearvars; tic;

%%% Choose one input method below: 
%%% Method 1 uses the UI, which pops up an window for user input.
%%% Method 2 is for setting the folders by command lines.
Method_of_choice =     2    ;

%%%%%======= Critical parameters for smFISH detection ==========
%%% Set the threshold for dot detection, the preset is 0.5
%%% Input values range [0 to infinity] but normally works 0.5 - 1.5.
thresForExon =      1.4     ; % <-- threshold for exon channels
thresForIntron =    1     ; % <-- threshold for intron channels
thresForNuc =       1     ; % <-- threshold for nuclear channels   

%%% indicate the type and order of the RNA channel to be analyzed.
%%% Put '1' for nucleus, '2' for exon, '3' for intron, '4' for protein IF,
%%% and '5' for skipping that channel.
%%% For example, if you have three channels in the order of exon Ch1, intron,
%%% and exon Ch2 then use: ChOrder = [ 2 1 2 ].
ChOrder = [  3  2  1  ];

%%% The parameters below are optimized for the C. elegans germline.
radius =     2.5     ; % the Radius of ROI for counting # dots per cell (microns)
vLim =       3       ; % maxium cell diameter for Voronoi cell segmentation.
nrange = [ 2.5  4.0 ]; % define range for diameter of nucleus (um).
sensi =      0.96    ; % sensitivity (circularity) for nucleus detection. The lower, the more stringent.
%%%%%============================================================= 

if Method_of_choice == 1
%%%%%---------- (Method 1) User Input (UI) system --------------------
    Image_path = uigetdir('C:\', 'Select a folder containing image files to process.');
    Image_path = strcat(Image_path, '\');
    Output_path = uigetdir(Image_path, 'Select a folder to save analysis results');
%%%%%--------------------------------------------------------------

elseif Method_of_choice == 2
%%%%%---------- (Method 2) Manual input for image folder location ---------
    Image_path1 = 'C:\CHL\smFISH images'; % <-- change folder as desired.
    Image_path2 = '';

    if Image_path1(end) ~= '\'
        Image_path1 = strcat(Image_path1, '\');
    elseif strcmp(Image_path1(end-1:end), '\\')
        Image_path1 = Image_path1(1:end-1);
    end
    Image_path = strcat(Image_path1, Image_path2);
    loc1 = strfind(Image_path, '\');
    dirName = Image_path(loc1(end-1)+1:end-1);
    

    Output_path1 = 'C:\CHL\smFISH analysis\temp'; % <-- change folder as desired.
    Output_path2 = '';  % <-- change folder as desired.

    if isempty(Output_path2)
        Output_path2 = dirName;
    end
    
    if Output_path1(end) ~= '\'
        Output_path1 = strcat(Output_path1, '\');
    elseif strcmp(Output_path1(end-1:end), '\\')
        Output_path1 = Output_path1(1:end-1);
    end
    
    if Output_path2(end) ~= '\'
        Output_path2 = strcat(Output_path2, '\');
    elseif strcmp(Output_path2(end-1:end), '\\')
        Output_path2 = Output_path2(1:end-1);
    end
    
    Output_path = strcat(Output_path1, Output_path2);
    loc2 = strfind(Output_path, '\');
    
    %%% Create an output folder to save results ---------
    for i = 1:length(loc2)
        if ~exist(Output_path(1:loc2(i)), 'dir')
            mkdir(Output_path(1:loc2(i)));
        end
    end
%%%%%---------------------------------------------------------------------

else
    fprintf('\nYou did not select the input method method (1 or 2)'\n);
end
savepath = strcat(Output_path, '\');
fprintf('\n\t\tThe radius of ROI is %2.1f um.\n\n', radius); 


%% Read and process image files
Imglist = dir(strcat(Image_path,'*.tif'));
tifVeri = 1;
if isempty(Imglist)
    Imglist = dir(strcat(Image_path,'*.lif'));
    tifVeri = 2;
end
if isempty(Imglist)
    tifVeri = 0;
    fprintf('\nThere is no image files in the folder.\n');
    return
end
Imglist = Imglist(~[Imglist.isdir]);
Imglist = {Imglist.name}';

af = cell(999,16);

%%% record the total # of images in the lif file
if tifVeri == 1
    fToRead = strcat(Image_path, Imglist);
    numOfImg = size(Imglist,1);
elseif tifVeri == 2 
    fToRead = strcat(Image_path, Imglist, '.lif');
    lifRead = bfGetReader(fToRead);
    omeMeta = lifRead.getMetadataStore();
    numOfImg = omeMeta.getImageCount();
end

%% ============= Read and process individual images ===============
for f=1:numOfImg     % open f-th image-stack in the lif file.
    %%% Use 'continue' if you want to exclude certain images from detection.        
    
    %         if f == 1 || f == 2
    %             continue
    %         end
    
    if tifVeri == 1
        imgdat = bfopen(fToRead{f});
    elseif tifVeri == 2
        imgdat = bfopenOne(fToRead,f);
    end

    %%% read metadata
    meta = imgdat{1, 4};
    Nch = meta.getPixelsSizeC(0).getValue();
    if Nch ~= length(ChOrder)
        Nch = length(ChOrder);
    end
    voxelSizeXdefaultValue = meta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
    voxelSizeXdefaultUnit = meta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
    voxelSizeX = meta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); 
    voxelSizeY = meta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); 
    voxelSizeZ = meta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); 

    iminfo = [
    Nch % # channel
    meta.getPixelsSizeX(0).getValue(); % image width, pixels
    meta.getPixelsSizeY(0).getValue(); % image height, pixels
    meta.getPixelsSizeZ(0).getValue(); % number of Z slices
    voxelSizeX.doubleValue();  % in µm
    voxelSizeY.doubleValue();  % in µm
    voxelSizeZ.doubleValue()   % in µm
    ];

    %% Separate channels
    %%% 'nucCh' contains all z-slices from the nucleus channel (ChOrder == 1)
    %%% 'exonCh1', 'exonCh2' from the exon channel (ChOrder == 2)
    %%% 'intronCh1', 'intronCh2' from the intron channel (ChOrder == 3)
    %%% 'elseCh' from the channel other than RNA or nuclei (ChOrder == 4)
    alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
    
    loc = 1:iminfo(1):iminfo(4)*iminfo(1);
    sigCount = 1;
    for i = 1:length(ChOrder)
        temp = imgdat{1,1}(loc+i-1,:);
        if ChOrder(i) == 1
            nucCh = temp;
        elseif ChOrder(i) == 2 || ChOrder(i) == 3 || ChOrder(i) == 4
            sigCh.(alph(sigCount)) = temp;
            sigCount = sigCount + 1;
        end
    end
    
    %%% Indicates how many exon/intron/IF channels are
    SigChOrder = ChOrder(ChOrder == 2 | ChOrder == 3 | ChOrder == 4);
    ChIndic = zeros(size(SigChOrder));
    for i = 2:4
        tempLoc = find(SigChOrder == i);
        ChIndic(tempLoc) = 1:length(tempLoc);
    end

    %% --------- Detect the gonad outline using an RNA channel --------
    if ismember(3, ChOrder)
        ChFor = sigCh.(alph(SigChOrder == 3 & ChIndic == 1));
    elseif ismember(2, ChOrder)
        ChFor = sigCh.(alph(SigChOrder == 2 & ChIndic == 1));
    else
        ChFor = nucCh;  % In case of no RNA channels, it uses the nuc channel.
    end

    pix = round(0.22 / mean(iminfo(5:6,1))); 
    if pix == 0
        pix = 1;
    end

    zproj = cat(3, ChFor{:,1});
    mip = max(zproj, [], 3);

    %%% "cyl_axis" contains: (unit: pixels)
    %%% (1) x-coor | (2) y-coor | (3) z-coor | (4) radius of the cell on that point.| 
    [cyl_axis, gonadBound, gonadMask] = DefOutlineATS(mip, pix, iminfo);

    %%% align germline orientation (distal end to the left)
    if mean(cyl_axis(1:10,4)) > mean(cyl_axis(end-9:end,4))
        cyl_axis(:,1) = iminfo(2)-cyl_axis(:,1)+1;
        gonadMask = flip(gonadMask, 2);
        imgdat{1,1}(:,1) = cellfun(@(x) flip(x,2), imgdat{1,1}(:,1), 'uni', 0);
    end

    
    %% --------------- Detect Nuclei -----------------------------
    %%% Detect nuclei and remove one nucleus of two overlapped nuclei.
    %%% NDNuc saves not overlapped Nuclei that are in DetNuc (detected Nuc).
    tnran = round(nrange/2 ./ iminfo(5:6,1)'); 
    tol = round( 1.2 / mean(iminfo(5:6,1)));  % tolerance for overlapped circles (um).
    
    %%% Set max pixel intensity value
    TifType = meta.getPixelsType(0).getValue();
    if strcmp(TifType, 'uint16')
        maxInt = 2^16 - 1;
    elseif strcmp(TifType, 'uint12')
        maxInt = 2^12 - 1;
    elseif strcmp(TifType, 'uint8')
        maxInt = 2^8 - 1;
    else
        maxInt = 2^8 - 1;
    end
        
    %%% Detect nuclei from image slices and make nuclear sphere. 
    %%% 'nucs' labels : | col 1: x-coor | y | z | radius | col 5: ave. circularity |
    %%%                 | col 6: std. circularity | col 7: total DAPI sig. area (# pixel) | 
    [Nuc, thrNuc, aftNuc] = DetectNucleus(nucCh, iminfo, pix, tnran, sensi, tol, f, numOfImg, thresForNuc, maxInt);

    %%% Quick visual
%     for i = 1:size(thrNuc,1)
%         figure,imshow(thrNuc{i,2})
%         pause
%         close all
%     end
    %%%--------------


    %%% Detect nuclei using bwconncomp and cross-validate the result with 
    %%% nuclei detected by hough transformation.
    %%% 'nucs' labels : | col 1: x-coor | y | z | radius | col 5: ave. circularity |
    %%%                 | col 6: std. circularity | col 7: total DAPI sig. area (# pixel) | 
    [nucs, DTCnuc] = DetectNucBlobs(thrNuc, aftNuc, Nuc, iminfo, tnran);

    %%% Remove nuclei that are outside the gonad
    nucs(:,8) = 1;
    for i = 1:size(nucs,1)
        r1 = round(nucs(i,2))-3;
        r2 = round(nucs(i,2))+3;
        r3 = round(nucs(i,1))-3;
        r4 = round(nucs(i,1))+3;
        if r1 < 1
            r1 = 1;
        end
        if r2 > iminfo(3)
            r2 = iminfo(3);
        end
        if r3 < 1
            r3 = 1;
        end
        if r4 > iminfo(2)
            r4 = iminfo(2);
        end
        %%% Use the flat gonad mask as a boundary for determining nucleus position relative to the gonad.
        if sum(sum(gonadMask(r1:r2, r3:r4))) == 0
            nucs(i,8) = 0;
        end
    end
    nucs(nucs(:,8) == 0,:) = [];
    nucs(:,8) = [];
    
    zproj = cat(3, nucCh{:,1});
    mip_nuc = max(zproj, [], 3);
    
    %%% Visualize detected nuclei
    %%%----------- max projection  ------------------
%         figure, imshow(mip_nuc*30);
%         hold on
%         plot(nucs(:,1), nucs(:,2), 'r.')
%         plot(Nuc(:,1), Nuc(:,2), 'ro')
    %%%----------------------------------------------


    %% ------------- Detect RNA spots ----------------------------
    %%% Store z-projected images of each signal channel in 'af'
    %%% 'loc': n-th column in 'af' to store the images 
    %%% (e.g. loc=10 to use column 10 and 11 to save images from 2 channels.
    loc = 10;  
    for i = 1:length(SigChOrder)
        zproj = cat(3, sigCh.(alph(i)){:,1});
        mip_Ch = max(zproj, [], 3);
        if i < 3
            af{f,loc+i-1} = mip_Ch;
        else
            loc = 10;
            af{f,loc+i-3} = mip_Ch;
        end
    end
    
    %%% process channels seperately for detecting RNA
    % 'blobrna': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
    %            | total intensity | x centroid | y centroid | z centroid.

    % 'blrna': | Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
    %          | 6th: matching nuc #(row # in 'nuc') 

    % 'rna' for intron Ch: | 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
    %               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
    %               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
    %               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

    % 'nuc for ATS' : | x-coor | y | z | radius | ave. circularity |
    %                 | 6th: std. circularity | total DAPI sig. | # trx sites | RNA ID .

    % 'nuc for mRNA' : | 1:3 xyz-coor | 4: radius | 5: ave. circularity |
    %           | 6: std. circularity | 7: total DAPI sig. area (# pixel) | 8: nuc ID.
    %           | 9: Voronoi mRNA count | 10: Voronoi count w/ limit (vLim, 3um) | 11: # mRNA in ROI w/ radius (2.5um)  |
    %  (aft_analysis->)       | 12th col: dist betw. nuc & cell outline | 13: 0 = cells on cortex, 1 = cells in rachis |

    %%% 'ExATSrna' column labels:
    %%% | 1:3 xyz-coor | 4: nuc ID containing the ATS | 5: ratio spot vol:mean single spot vol | 
    %%% | 6: raw intensity from exon Ch | 7: center z-plane | 8: RNA ID | 
    %%% | 9: # mRNA equivalent (by sig. intensity) |
    
    %%% 'ExATSnuc'
    %%% 'nuc for ATS' : | 1-3: xyz-coor | 4: radius (um) | 5: ave. circularity |
    %%%                 | 6: std. circularity | 7: total DAPI sig. | 8: # ATS | 
    %%%                 | 9: total raw ATS intensity (all ATSs) | 10: total # mRNA equivalent per nuc|

    pixr = round( 0.1 / mean(iminfo(5:6,1)));
    rna = cell(length(SigChOrder),1);
    nuc = cell(length(SigChOrder),1);
    
    parfor g = 1:length(SigChOrder)
        alphPar = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
        gmask = gonadMask;
        rZproj = af{f,9+g};
        rPla = sigCh.(alphPar(g))(:,1);
        %%%%%------------- Process EXON channels -----------------------
        if SigChOrder(g) == 2   % 2 is for processing exon images.
            [DETrna, derna, mom2(g)] = DetectRNAExon(rPla, pixr, iminfo, f, numOfImg, g, thresForExon);
            if ~isempty(DETrna)
                %%% Put results of RNA1 and RNA2 in separate struct.
                DETrna(:,8) = 1:length(DETrna(:,1)); % mRNA ID.
                [mrnas, mrnasNoATS, rnas_atsEx, nucm, nuc_atsEx] = MatchmrnaNuc(DETrna, thrNuc, nucs, iminfo, radius, vLim, pixr);
                
                %%% Remove RNA outside the gonad
                mrnas(:,13) = 0;
                for i = 1:size(mrnas,1)
                    if gmask(round(mrnas(i,2)), round(mrnas(i,1))) == 0
                        mrnas(i,13) = 999;
                    end
                end
                mrnas(mrnas(:,13) == 999,:) = [];
                mrnas(:,13) = [];
                
                rnas_atsEx(:,10) = 0;
                for i = 1:size(rnas_atsEx,1)
                    if gmask(round(rnas_atsEx(i,2)), round(rnas_atsEx(i,1))) == 0
                        rnas_atsEx(i,13) = 999;
                    end
                end
                rnas_atsEx(rnas_atsEx(:,10) == 999,:) = [];
                rnas_atsEx(:,10) = [];
                
                if ismember(3, ChOrder)
                    rna{g} = mrnas; % save info for detected mRNA
%                     rna{g}(:,4) = 1;  % Assuming all mRNA are paired to nucleus.
                else
                    rna{g} = mrnasNoATS;
                end
                nuc{g} = nucm;

                %%% record detected ATS from exon channel
                ExATSrna{g} = rnas_atsEx;
                ExATSnuc{g} = nuc_atsEx;
            else
                rna{g} = [];
                nuc{g} = [];
            end
            
            %%%---- Quick visual ----
%             figure,imshow(af{f,9+g}*100)
%             figure,imshow(af{f,9+g}*100)
%             hold on
%             plot(rna{g}(:,1), rna{g}(:,2), 'ro');
            %%%-----------------------

        %%%%%-------------- Process INTRON channels --------------------   
        elseif SigChOrder(g) == 3    % 3 is for processing intron images.
            % mom1 or mom2 is the mean background signal intensity from all z-planes in a gonad.
            [blobrna, blrna, mom1(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, numOfImg, g, thresForIntron);

            if isempty(blobrna) == 0
                % determine which nuclei RNA spots belong. 
                [prnas, nuct] = MatchrnaNuc(thrNuc, blobrna, nucs, iminfo, pixr, 10);
                
                %%% Remove dim RNA spots using signal-to-local background ratio on z-projected image.
                prnas(:,8) = 0;
                ratioSB = zeros(size(prnas,1),1);
                if ~isempty(prnas)
                    for j = 1:size(prnas,1)
                        rng1 = pixr*2;
                        rng2 = pixr*10;
                        
                        Xrange1 = round(prnas(j,1))-rng1:round(prnas(j,1))+rng1;
                        Yrange1 = round(prnas(j,2))-rng1:round(prnas(j,2))+rng1;
                        Xrange1 = Xrange1(Xrange1 > 0 & Xrange1 <= iminfo(2));
                        Yrange1 = Yrange1(Yrange1 > 0 & Yrange1 <= iminfo(3));
                        sigMean1 = rZproj( Yrange1 , Xrange1);
                        sigMean2 = mean(mean(sigMean1));
                        
                        Xrange2 = ceil(prnas(j,1))-rng2:ceil(prnas(j,1))+rng2;
                        Yrange2 = ceil(prnas(j,2))-rng2:ceil(prnas(j,2))+rng2;
                        Xrange2 = Xrange2(Xrange2 > 0 & Xrange2 <= iminfo(2));
                        Yrange2 = Yrange2(Yrange2 > 0 & Yrange2 <= iminfo(3));
                        bgMean1 = rZproj( Yrange2 , Xrange2 );
                        bgMean1(rng2-rng1*2+1:rng1+rng2+1, rng2-rng1*2+1:rng1+rng2+1) = 0;
                        bgMean2 = bgMean1(bgMean1 > 0);
                        sortN2 = sort(bgMean2(:));
                        bgMean3 = mean(sortN2(10:end));
                        ratioSB(j) = sigMean2/bgMean3;
%                         fprintf('\n%.2f for %d', ratioSB(j), j)
                        if ratioSB(j) < 1.05   % the cutoff signal-local bg ratio of each mRNA
                            prnas(j,8) = 999;
                        end
                    end
                end
                prnas(prnas(:,8) == 999,:) = [];
                prnas(:,8) = [];

%%%----------------------- Quick visual ---------------------------
%                 figure, imshow(rZproj*50)
%                 hold on
% %                 plot(blobrna(:,7), blobrna(:,8), 'r+');
%                 plot(prnas(:,1), prnas(:,2), 'ro');
%                 for i = 1:size(blobrna,1)
% %                     text(blobrna(i,7), blobrna(i,8), num2str(i), 'color', 'c');
%                     text(prnas(i,1), prnas(i,2), num2str(i), 'color', 'c');
%                 end
%%%-----------------------------------------------------------------
                
                % Put results of RNA1 and RNA2 in separate struct.
                rna{g} = prnas;
                nuc{g} = nuct;
            else
                rna{g} = [];
                nuc{g} = [];
            end

        %%%%%-------------- Process other channels --------------------
        elseif SigChOrder(g) == 4   
            %%% Implement for IF analysis
        end
    end 

%%%----------- Quick visual for RNA detection validation --------
%     ch = 1; % The channel to visualize
%     figure,imshow(af{f,ch+9}*100)
%     figure,imshow(af{f,ch+9}*100)
%     hold on
%     plot(rna{ch}(:,1), rna{ch}(:,2), 'r.')
%%%---------------------------------------------------------------


    %% Cross-validate ATS between intron and exon channels.
    % Also remove mRNA spots that are determined as ATS.

    % 'DexRext' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
    %               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
    %               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
    %               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

    if ismember(2, SigChOrder)
        if ~isempty(rna{SigChOrder == 2})
            mRNAIntM = mean(rna{SigChOrder == 2}(:,6)); % mean intensity of all mRNA spots.
        else
            mRNAIntM = 0.0001;
        end
    else
        mRNAIntM = 0.0001;
    end

    if ismember(3, SigChOrder)
        if ~isempty(rna{SigChOrder == 3})
            pRNAIntM = mean(rna{SigChOrder == 3}(:,6)); 
        else
            pRNAIntM = 0.0001;
        end
    else
        pRNAIntM = 0.0001;
    end

    if ismember(2, SigChOrder) && ismember(3, SigChOrder)
        rna{SigChOrder == 2}(:,13) = 0;
        rna{SigChOrder == 3}(:,8:9) = 0;
        
        distATS =  1.5  ;   % acceptible 3D distance (um) between ATSs from the two channels.

        countp = 0;
        countm = 0;

        % remove detected ATS spots from the list of mRNAs
        if ~isempty(rna{SigChOrder == 3})
            for i = 1:length(rna{SigChOrder == 3}(:,1))
                zpln = rna{SigChOrder == 3}(i,7);
                DetRext = rna{SigChOrder == 2}(rna{SigChOrder == 2}(:,7) > zpln- 6 & rna{SigChOrder == 2}(:,7) < zpln+ 6,:);
                distMat = zeros(length(DetRext(:,1)),5);
                distMat(:,1) = DetRext(:,1) - rna{SigChOrder == 3}(i,1);   
                distMat(:,2) = DetRext(:,2) - rna{SigChOrder == 3}(i,2);  
                distMat(:,3) = DetRext(:,6);
                distMat(:,4) = sqrt(distMat(:,1).^2 + distMat(:,2).^2); % + distMat(:,3).^2);   % X,Y,Z distances are all used.
                distMat(:,5) = DetRext(:,8);
                distMats = distMat(distMat(:,4) <= distATS/iminfo(6), :);
                if isempty(distMats)
                    rna{SigChOrder == 3}(i,9) = 999;
                    continue
                end
                distMats = sortrows(distMats,-3);

                if distMats(1,4) <= distATS/iminfo(6)
                    if distMats(1,3) < mRNAIntM/2 && rna{SigChOrder == 3}(i,6) < pRNAIntM/5 && rna{SigChOrder == 3}(i,1) > min(rna{SigChOrder == 3}(:,1)) + 25/iminfo(6)
                        rna{SigChOrder == 3}(i,9) = 999;  % mark intron spot unlikely ATS.
                        countp = countp + 1;
                    else
                        rna{SigChOrder == 2}(rna{SigChOrder == 2}(:,8) == distMats(1,5),13) = 999;  % mark ATS in exon channel
                        rna{SigChOrder == 3}(i,8) = distMats(1,3);
                        countm = countm + 1;
                    end
                else
                    rna{SigChOrder == 3}(i,9) = 999;
                end
            end
        end

        % remove intron ATS that are not seen at exon channel
%         NoOvlapATS = prna{ChOrder == 3}(prna{ChOrder == 3}(:,8) == 999,:);
        if countp > 0
%                 prna{ChOrder == 3}(prna{ChOrder == 3}(:,9) == 999,:) = [];
        end
        rna{SigChOrder == 3}(:,9) = [];

        if countm > 0
%                 prna{ChOrder == 2}(prna{ChOrder == 2}(:,13) == 999,:) = [];  
        end
        rna{SigChOrder == 2}(:,13) = [];
    end



    %%%%%---------- Record ATS intensity from exon channel ------------
    %%% 'nuc' for ATS
    %       |1-3: XYZ-coordinates |4: radius (pixels)| 5: ave. circularity| 6: std circularity |
    %       | 7: DAPI intensity |8: # ATS per nuc | 9: total ATS intensity from INTRON channel |
    %       | 10: total ATS int. from EXON channel |
    %       | 11: # mRNA equivalent for ATS |
    if ismember(2, SigChOrder) && ismember(3, SigChOrder)
        nNcount = find( nuc{SigChOrder == 3}(:,8) > 0 );
        for i = 1:length(nNcount)
            ints = rna{SigChOrder == 3}(rna{SigChOrder == 3}(:,4) == nNcount(i),6);
            nuc{SigChOrder == 3}(nNcount(i),9) = sum(ints);
            ints2 = rna{SigChOrder == 3}(rna{SigChOrder == 3}(:,4) == nNcount(i),8);
            nuc{SigChOrder == 3}(nNcount(i),10) = sum(ints2);
            nuc{SigChOrder == 3}(nNcount(i),11) = nuc{SigChOrder == 3}(nNcount(i),10)/mRNAIntM;
        end
    end

    %% Make summary data
    % 'af': | 1: nuc info (Ch1) | 2: RNA info (Ch1) | 3: nuc info (Ch2) | 4: RNA info (Ch2) | .
    %       | 5: cyl_axis info | 6: gonad outline mask | 7: z-projection of nucleus channel |.
    %       | 8: mean background intensity from Ch1 | 9: mean bg. noise from Ch2 |
    %       | 10: z-proj of Ch1 | 11: z-proj of Ch2 |
    %       (after analysis)     
    %       | 12th col: distance distal tip & first chain nuc | 12: # nuc in rachis |

    % record DTC? Use locator.
    
    afLoc = 1;
    for i = 1:length(SigChOrder)
        if ~isempty(rna{i})
            temp = rna{i}(rna{i}(:,4) ~= 0,:);
            if i < 3
                loc = i;
                af(f,loc*2-1) = nuc(i);
                af(f,loc*2) = {temp};
            else
                loc = 12;   % To put Ch3 info in col 17 & 18 
                af(f,i*2-1 + loc) = nuc(i);
                af(f,i*2 + loc) = {temp};
            end
        end
    end

    af{f,5} = cyl_axis;
    af{f,6} = gonadMask;
    af{f,7} = mip_nuc;
    if exist('mom1', 'var')
        af{f,8} = mom1(mom1~=0);
    end
    if exist('mom2', 'var')
        af{f,9} = mom2(mom2~=0);
    end

    %%% af: col 10 and 11 are for z-proj images of Ch1 and Ch2.
    %%% af: | col 12: nuc info | 13: ATS info | <-- using only exon channel
    %%%     | 14: nuc info (exon Ch2) | 15: ATS info (exon Ch2) if two exon channels exist.
    if exist('ExATSrna', 'var') && ismember(2, SigChOrder)
        loc = 12;
        for i = 1:length(ExATSrna)
            if ~isempty(ExATSrna{i})
                temp = ExATSrna{i}(ExATSrna{i}(:,4) ~= 0,:);
                af(f,loc) = ExATSnuc(i);
                af(f,loc+1) = {temp};
                loc = loc + 2;
            end
        end
    end
    
    %%% af: | col 16: 'iminfo' | 17: info for DTC nucleus |
    af{f,16} = iminfo;
    af{f,17} = DTCnuc;
    
    %%% if 3rd RNA/IF channel exists
    %%% af: | 18: nuc info (Ch3) | 19: RNA info (Ch3) |

%     cd(savepath);
%     save(strcat('smFISHworkspace_middle_',num2str(f)), '-regexp', '^(?!(imgdat|lifRead|meta|omeMeta)$).');

end

af(numOfImg+1:end,:) = [];

cd(savepath);
save('smFISHworkspace_final', '-regexp', '^(?!(imgdat|lifRead|meta|omeMeta)$).');

eTime = toc;
fprintf('\n\t\tElapsed time is %2.1f minutes.\n\n', eTime/60); 

