classdef Prj2TB
    % Functions used for ECES 435 Project 2
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Prj2TB(self)
            % imgs = Prj2TB.read_all_jpgs('raw_img/FV5', '', 'TIF');
            % imgs_cropped = Prj2TB.crop_images(imgs);
            
        end
    end
    
    methods (Static)
        
        function img_out = read_all_imgs(dir_path, prefix, ext)
            dir_s = dir(dir_path);
            img_out = {};
            regexp_pat = sprintf('%s\.%s$',prefix,ext);
            for i = dir_s'
                hit = regexpi(i.name, regexp_pat);
                if ~isempty(hit)
                    try
                        path = fullfile(dir_path, i.name);
                        img_out{end+1} = imread(path);
                        fprintf('loaded: %s\n', path);
                    catch
                        fprintf('*** ERROR loading: %s\n', path);
                    end
                end
            end            
        end
        
        function img_out = crop_images(img_in) 
            img_out = {};
            for i = img_in
                img = i{1};
                img_size = size(img);
                center = floor(img_size(1:2)/2);
                box = floor([1024 1024]/2)*2;
                lc = center - box/2 + 1;
                uc = center + box/2;
                img_out{end+1} = img(lc(1):uc(1), lc(2):uc(2), :);            
            end            
        end
        
        function [] = save_images(img_in, dir, prefix)
            len = length(img_in)
            for i = 1:len
                filename = sprintf('%s_%03i.TIF', prefix, i);
                filepath = fullfile(dir, filename)
                imwrite(img_in{i}, filepath);
            end
        end
        
        function img_aligned = align(imgA_raw, imgB_raw, sw_plot)
            ptThresh = 0.15;
            
            if sw_plot
                figure; 
                imshowpair(imgA_raw, imgB_raw, 'montage');
                title('original images');
            end
            
            SS = 1;
            % select Y channel
            imgA = imresize(rgb2ycbcr(imgA_raw), SS);
            imgB = imresize(rgb2ycbcr(imgB_raw), SS);
            imgA = imgA(:,:,1);
            imgB = imgB(:,:,1);
            
            pointsA = detectFASTFeatures(imgA, 'MinContrast', ptThresh/SS^.5);
            pointsB = detectFASTFeatures(imgB, 'MinContrast', ptThresh/SS^.5);
            
            %{
            if sw_plot
                figure;
                subplot(1,2,1); imshow(imgA); hold on;
                plot(pointsA); hold off;
                title('Corners in A');
                subplot(1,2,2); imshow(imgB); hold on;
                plot(pointsB); hold off;
                title('Corners in B');
            end
            %}
            
            [featuresA, pointsA] = extractFeatures(imgA, pointsA, 'BlockSize', 1+10*SS);
            [featuresB, pointsB] = extractFeatures(imgB, pointsB, 'BlockSize', 1+10*SS);
            indexPairs = matchFeatures(featuresA, featuresB, 'MatchThreshold', SS*10);
            pointsA_matched = pointsA(indexPairs(:, 1), :);
            pointsB_matched = pointsB(indexPairs(:, 2), :);
            
            shift = pointsB_matched.Location - pointsA_matched.Location;
            
            % prepare for DBSCAN
            % find 2 norm (distance to origin)
            % dto = sqrt( sum( (shift.^2)' )' );
            % epsilon = median(dto);
            epsilon = 2 * SS^2;
            min_pts = ceil(size(shift, 1) * 0.05);
            idx = DBSCAN(shift, epsilon, min_pts);
            shift_sel = shift(idx==1, :);
            pointsA_sel = pointsA_matched(idx == 1);
            pointsB_sel = pointsB_matched(idx == 1);
            
            fprintf('num of corners matched: %d\n',length(shift));
            fprintf('num of corners selected: %d\n',sum(idx==1));
            
            
            
            if sw_plot
                figure;
                subplot(1,2,1); imshow(imgA); hold on;
                plot(pointsA); plot(pointsA_matched.Location(:,1),pointsA_matched.Location(:,2), 'o'); 
                hold off; title('Corners in A');
                subplot(1,2,2); imshow(imgB); hold on;
                plot(pointsB); plot(pointsB_matched.Location(:,1),pointsB_matched.Location(:,2), 'o'); 
                hold off; title('Corners in B');
                
                
                figure; 
                showMatchedFeatures(imgA, imgB, pointsA_matched, pointsB_matched);
                legend('A', 'B');
                title('matched features comparison (overlay)')
                
                figure; 
                showMatchedFeatures(imgA, imgB, pointsA_matched, pointsB_matched, 'montage');
                legend('A', 'B');
                title('matched features comparison (side-by-side)')
                
                figure; hold on;
                plot(shift(:,1),shift(:,2), '.')
                plot(shift_sel(:,1),shift_sel(:,2), 'o')
                grid on; hold off; axis square
                title('scatter plot of feature shifts for clustering')
                legend('shift vecotr','shift vecotr selected');
            end
            
            [tform, pointsBm, pointsAm] = estimateGeometricTransform(...
                pointsB_sel, pointsA_sel, 'affine');
            img_aligned = imwarp(imresize(imgB_raw,SS), tform, 'OutputView', imref2d(size(imgB)));
            % pointsBmp = transformPointsForward(tform, pointsBm.Location);
            
        end
            
        
    end
        
end

