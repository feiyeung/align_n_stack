%{
%% load and crop
imgs = Prj2TB.read_all_imgs('raw_img/FV5','','TIF');
imgs_cropped = Prj2TB.crop_images(imgs);

%% show
close all
pause(0.5); figure;
for i=1:length(imgs_cropped)
    imshow(imgs_cropped{i});
    title(num2str(i));
    pause(3);
end
%}





%% clear and load cropped images
clear;
imgs = Prj2TB.read_all_imgs('.','','TIF');
close all


%% align all?
close all
imgs_tfed = {};
imgs_cumu = {};
ind_range = [1,2,5,6,7,8,9];
ind_range = [5,6,7,8,9];
% ind_range = [5,6];
imtermediate_plots = false;
count = 1;
img_base = imgs{ind_range(1)};
img_pool = uint32(imresize(img_base,1));
for i = ind_range(2:end)
    count = count + 1;
    if count>3 && imtermediate_plots
        close all
    end
    this_img = Prj2TB.align(img_base, imgs{i}, imtermediate_plots);
    size(this_img)
    size(img_pool)
    imgs_tfed{end+1} = this_img;
    img_pool = img_pool + uint32(this_img);
    imgs_cumu{end+1} = uint16(img_pool / count);   
    figure;
    imshow(imgs_cumu{end});
    title(sprintf('img number %i', i))
    
    %{
    center = flip([509, 780]);
    box = floor([256 256]/2)*2;
    lc = center - box/2 + 1;
    uc = center + box/2;
    this_img = imgs_cumu{end}(lc(1):uc(1), lc(2):uc(2), :);       
    filename = sprintf('progressive_%i.PNG',count);
    imwrite(imresize(uint8(this_img ./ 2^8), 8, 'nearest'), filename);
    %}
end
img_out = uint16(img_pool ./ length(ind_range));

%{
this_img = img_base(lc(1):uc(1), lc(2):uc(2), :);       
filename = 'progressive_0.PNG';
imwrite(imresize(uint8(this_img ./ 2^8), 8, 'nearest'), filename);

%}

%% show and compare

figure; imshowpair(img_out, img_base, 'montage');




figure; imshowpair(img_out(lc(1):uc(1), lc(2):uc(2), :), img_base(lc(1):uc(1), lc(2):uc(2), :), 'montage');


%%
close all

imgA_raw = imgs{1};
imgB_raw = imgs{2};

pause(0.5); figure; imshow(imgA_raw);
title('original image A');

pause(0.5); figure; imshow(imgA_raw);
title('original image B');


imgA = imgA_raw(:,:,1);
imgB = imgB_raw(:,:,1);

ptThresh = 0.2;
pointsA = detectFASTFeatures(imgA, 'MinContrast', ptThresh);
pointsB = detectFASTFeatures(imgB, 'MinContrast', ptThresh);

pause(0.5); figure; imshow(imgA); hold on;
plot(pointsA);
title('Corners in A, red channel');

pause(0.5); figure; imshow(imgB); hold on;
plot(pointsB);
title('Corners in B, red channel');

%%
[featuresA, pointsA] = extractFeatures(imgA, pointsA);
[featuresB, pointsB] = extractFeatures(imgB, pointsB);

%%
indexPairs = matchFeatures(featuresA, featuresB);
pointsA = pointsA(indexPairs(:, 1), :);
pointsB = pointsB(indexPairs(:, 2), :);



%%
pause(0.5); figure; showMatchedFeatures(imgA, imgB, pointsA, pointsB);
legend('A', 'B');
title('features marked and paired between A and B')

%%
shift = pointsA.Location - pointsB.Location;
pause(0.5); figure;
plot(shift(:,1),shift(:,2), 'x')
grid on;
axis([-200 200 -200 200])
title('scatter plot of feature shifts for clustering')
axis square

pt_ref = [-16 * ones(length(shift), 1),-23 * ones(length(shift), 1)] ;
norm2sq = sum(transpose((shift - pt_ref).^2));
sel_ind = norm2sq < 20^2;

ptA_sel = pointsA(sel_ind);
ptB_sel = pointsB(sel_ind);


%% change to selected points
pointsA = ptA_sel;
pointsB = ptB_sel;


%%
[tform, pointsBm, pointsAm] = estimateGeometricTransform(...
    pointsB, pointsA, 'affine');
imgBp = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB)));
pointsBmp = transformPointsForward(tform, pointsBm.Location);

%%
pause(0.5); figure;
showMatchedFeatures(imgA, imgBp, pointsAm, pointsBmp);
legend('A', 'B');
title('image B shifted and overlayed with corners marked')


%% 

imgA = imgA_raw;
imgB = imgB_raw;

% Extract scale and rotation part sub-matrix.
H = tform.T;
R = H(1:2,1:2);
% Compute theta from mean of two possible arctangents
theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
% Compute scale from mean of two stable mean calculations
scale = mean(R([1 4])/cos(theta));
% Translation remains the same:
translation = H(3, 1:2);
% Reconstitute new s-R-t transform:
HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; ...
  translation], [0 0 1]'];
tformsRT = affine2d(HsRt);

imgBold = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB)));
imgBsRt = imwarp(imgB, tformsRT, 'OutputView', imref2d(size(imgB)));

pause(0.5); figure, clf;
imshowpair(imgBold,imgBsRt,'ColorChannels','red-cyan'), axis image;
title('Color composite of affine and s-R-t transform outputs');

%%
pause(0.5); figure;
showMatchedFeatures(imgA, imgBsRt, pointsAm, pointsBmp);
legend('A', 'B');
title('image B shifted and overlayed with corners marked')

%% 
close all
im_merge = (uint32(imgA) + uint32(imgBsRt))./2;
im_merge = uint16(im_merge);

pause(0.5); figure;
imshow(imgA)
title('original image A');

pause(0.5); figure;
imshow(im_merge)
title('shifted and stacked with mean value');

img_out = {};
for i = {im_merge, imgA}
    this_img = i{1};
    img_size = size(this_img)
    % center = floor(img_size(1:2)/2) 
    center = [180, 780];
    box = floor([256 256]/2)*2;
    lc = center - box/2
    uc = center + box/2
    img_out{end+1} = this_img(lc(1):uc(1), lc(2):uc(2), :);            
end      
pause(0.5); figure; imshowpair(img_out{1}, img_out{2}, 'montage');
title('256px cropped comparison bettwen stacked(L) and raw(R)');


