function outp = lifopen(path)

imdat = bfopen('D:\smFISH images\prac2.lif');
% subject = metadata.get('Subject');
% title = metadata.get('Title');



%% Retrieving metadata

% Original metadata is a set of key/value pairs specific to the input format of the data. It is stored in the data{s, 2} element of the data structure returned by bfopen.
metadata = lifdat{1, 2};
metadataKeys = metadata.keySet().iterator();
for i=1:metadata.size()
  key = metadataKeys.nextElement();
  value = metadata.get(key);
  fprintf('%s = %s\n', key, value)
end

% OME metadata is a standardized metadata structure, which is the same regardless of input file format. It is stored in the data{s, 4} element of the data structure returned by bfopen, and contains common metadata values such as physical pixel sizes, instrument settings, and much more. See the OME Model and Formats documentation for full details.
imn = 1; % image number
meta = imdat{imn, 4};
imn = imn - 1;
iminfo = [
meta.getPixelsSizeX(imn).getValue(); % image width, pixels
meta.getPixelsSizeY(imn).getValue(); % image height, pixels
meta.getPixelsSizeZ(imn).getValue(); % number of Z slices
meta.getPixelsPhysicalSizeX(imn).getValue(); % in ?m
meta.getPixelsPhysicalSizeY(imn).getValue(); % in ?m
meta.getPixelsPhysicalSizeZ(imn).getValue(); % in ?m
meta.getChannelCount(imn); % # channel
];
