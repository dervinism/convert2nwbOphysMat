function images = readMultiLayerTiff(imageFilename)

info = imfinfo(imageFilename);
numberOfImages = length(info);
for img = 1:numberOfImages
    images{img} = imread(imageFilename, img); %#ok<*AGROW> 
end	