set(0,'DefaultFigureWindowStyle','docked')

% load images
reporter = double(imread('p2_T052_C0.tif'));
DIC = double(imread('p2_T052_C1.tif'));
nuc = double(imread('p2_T052_C2.tif'));

% Display images
figure('name','DIC','NumberTitle', 'off'); imshow(DIC,[])
figure('name','nuc','NumberTitle', 'off'); imshow(nuc,[])
figure('name','reporter','NumberTitle', 'off'); imshow(reporter,[])

%% Create mask of cells
% smooth
nuc_smooth = imgaussfilt(nuc,2); % maybe don't smooth for the mask
figure('name','nuc_smooth','NumberTitle', 'off'); imshow(nuc_smooth,[])

% Finding the best threshold manually :-(, could use Otsu! :-)
% for x=20:5:40
%   nuc_thresh = nuc_smooth > x; % static threshold
%   figure('name',num2str(x),'NumberTitle', 'off'); imshow(nuc_thresh,[])
% end

% threshold
nuc_thresh = nuc_smooth > 20; % static threshold
%nuc_thresh = nuc_smooth > graythresh(nuc_smooth); % otsu dynamic threshold
figure('name','nuc_thresh','NumberTitle', 'off'); imshow(nuc_thresh,[])

% clear small objects
nuc_nuc_cleared = bwareaopen(nuc_thresh, 2000);
figure('name','nuc_nuc_cleared','NumberTitle', 'off'); imshow(nuc_nuc_cleared,[])

% fill holes
nuc_filled = imfill(nuc_nuc_cleared,'holes');
figure('name','nuc_filled','NumberTitle', 'off'); imshow(nuc_filled,[])

% rename variable
nuc_mask = nuc_nuc_cleared;

% % Finding the best number manually :-(, could use Otsu! :-)
% for x=10:3:30
%   nuc_smooth = imgaussfilt(nuc,x);
%   seeds = imregionalmax(nuc_smooth);
%   seeds(nuc_mask==0)=0;
%   [X Y] = find(seeds);
%   figure('name',num2str(x),'NumberTitle', 'off')
%   imshow(nuc,[]);
%   hold on;
%   plot(Y,X,'or','markersize',2,'markerfacecolor','r')
% end


%% Find seeds
nuc_smooth = imgaussfilt(nuc,19); % more smoothing to join close seeds
seeds = imregionalmax(nuc_smooth);
seeds(nuc_mask==0)=0; % remove seeds outside of our nuc mask
% Debug with plot
[X Y] = find(seeds);
figure('name','seeds','NumberTitle', 'off')
imshow(nuc,[]);
hold on;
plot(Y,X,'or','markersize',2,'markerfacecolor','r')


%% Watershed
nuc_smooth = imgaussfilt(nuc,1); % don't smooth too much for watersheding
nuc_min = imimposemin(-nuc_smooth,seeds); % set locations of seeds to be -Inf (because matlab's watershed uses that to mark seperate segments)
figure('name','nuc_min','NumberTitle', 'off'); imshow(nuc_min,[]);
nuc_ws = watershed(nuc_min);
figure('name','nuc_ws','NumberTitle', 'off'); imshow(nuc_ws,[]);
nuc_ws(nuc_mask==0)=0; % remove areas that aren't in our nuc mask
figure('name','nuc_ws','NumberTitle', 'off')
imshow(nuc_ws,[]);

% rename variable
labelled_nuc = nuc_ws; % set a better variable name

% Clear cells touching the boarder
labelled_nuc = imclearborder(labelled_nuc);

%% Color overlay v1
% Display original image
figure
imshow(imoverlay(uint8(nuc),bwperim(labelled_nuc),'r'),[]);
hold on
% Display color overlay
labelled_nuc_rgb = label2rgb(uint32(labelled_nuc), 'jet', [1 1 1], 'shuffle');
himage = imshow(labelled_nuc_rgb,[]);
himage.AlphaData = nuc_mask*.3;
% Display green seed dots
[xm,ym]=find(seeds);
hold on
plot(ym,xm,'or','markersize',2,'markerfacecolor','g','markeredgecolor','g')



%% Color overlay v2 (perim only)
% Display original image
figure
imshow(imoverlay(uint8(nuc),bwperim(labelled_nuc),'r'),[]);
hold on
% Display green seed dots
[xm,ym]=find(seeds);
hold on
plot(ym,xm,'or','markersize',2,'markerfacecolor','g','markeredgecolor','g')


%% Bar graphs per cell
% Calculate total intensity per cell
total_intensity = [];
for id=unique(labelled_nuc(:))'
    if id==0 % skip background
        continue
    end
    single_cell_mask= labelled_nuc==id;
    %imshow(single_cell_img,[])
    single_cell_values = nuc(single_cell_mask)
    total_intensity = [total_intensity sum(single_cell_values)] 
end

% Calculate area per cell
stats = regionprops(labelled_nuc,nuc,'Area');

% Plot total intesity bar graph
figure
bar(total_intensity)
title('total intensity')

% Plot area bar graph
area = cat(1,stats.Area);
figure
title('area')
bar(area)

% Plot area with total intensity
bar([area(2:13)'; total_intensity])


% Image processing experiment on DIC image
im=DIC(1:2333,1:2333); % crop
im=imgradient(im); % first derivative image
im_smooth = imgaussfilt(im,2);
figure('name','im','NumberTitle', 'off'); imshow(im_smooth,[])
im_thresh = im_smooth>1000;
figure('name','im_thresh','NumberTitle', 'off'); imshow(im_thresh,[])
% fill holes
nuc_filled = imfill(im_thresh,'holes');
figure('name','nuc_filled','NumberTitle', 'off'); imshow(nuc_filled,[])
% clear small objects
nuc_nuc_cleared = bwareaopen(nuc_filled, 2000);
figure('name','nuc_nuc_cleared','NumberTitle', 'off'); imshow(nuc_nuc_cleared,[])
