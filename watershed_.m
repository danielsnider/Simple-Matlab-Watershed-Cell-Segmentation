set(0,'DefaultFigureWindowStyle','docked')

cyto = imread('Z:\OPRETTA\Operetta Raw Data\Mammalian cells\20171020_LFS_DPC__2017-10-20T18_27_46-Measurement2\Images\r03c06f01p01-ch1sk6fk1fl1.tiff');
cyto = imread('\\carbon.research.sickkids.ca\rkafri\OPRETTA\Operetta Raw Data\Mammalian cells\20171020_LFS_DPC__2017-10-20T18_27_46-Measurement2\Images\r03c06f02p01-ch1sk37fk1fl1.tiff'); % dim image

figure;
imshow(cyto,[prctile(cyto(:),1) prctile(cyto(:),99)]);


%% Remove extreme pixel values (this also scales values between 0 and 1)
bot = double(prctile(cyto(:),1));
top = double(prctile(cyto(:),99));
cyto = (double(cyto)-bot) / (top - bot);
cyto(cyto>1) = 1; % do the removing
cyto(cyto<0) = 0; % do the removing
figure('name','cyto','NumberTitle', 'off'); imshow(cyto,[])

%% Create mask of cells
% smooth
cyto_smooth = imgaussfilt(cyto,2); % maybe don't smooth for the mask
figure('name','cyto_smooth','NumberTitle', 'off'); imshow(cyto_smooth,[])

% threshold
cyto_thresh = cyto > .18;
%cyto_thresh = cyto > graythresh(cyto); % otsu
figure('name','cyto_thresh','NumberTitle', 'off'); imshow(cyto_thresh,[])

% clear small objects
cyto_nuc_cleared = bwareaopen(cyto_thresh, 2000);
figure('name','nuc_mask_cleared','NumberTitle', 'off'); imshow(cyto_nuc_cleared,[])
% % Finding the best number manually :-(, could use Otsu! :-)
% for x=800:100:2000
%     nuc_mask_cleared = bwareaopen(cyto_thresh, x);
%     %nuc_mask = imopen(cyto_thresh,strel('disk',x)); 
%     figure('name',num2str(x),'NumberTitle', 'off'); imshow(nuc_mask_cleared,[])
% end

cyto_mask = cyto_nuc_cleared; % set a better variable name


%% Find seeds
cyto_smooth = imgaussfilt(cyto,12); % more smoothing to join close seeds
seeds = imregionalmax(cyto_smooth);
seeds(cyto_mask==0)=0; % remove seeds outside of our cyto mask
% Debug with plot
[X Y] = find(seeds)
figure('name','seeds','NumberTitle', 'off')
imshow(cyto,[]);
hold on;
plot(Y,X,'or','markersize',2,'markerfacecolor','r')

% Finding the best number manually :-(, could use Otsu! :-)
% for x=10:20
%   cyto_smooth = imgaussfilt(cyto,x);
%   seeds = imregionalmax(cyto_smooth);
%   seeds(cyto_mask==0)=0;
%   [X Y] = find(seeds)
%  figure('name',num2str(x),'NumberTitle', 'off')
%   imshow(cyto,[]);
%   hold on;
%   plot(Y,X,'or','markersize',2,'markerfacecolor','r')
% end


%% Watershed
cyto_smooth = imgaussfilt(cyto,1); % don't smooth too much for watersheding
cyto_min = imimposemin(-cyto_smooth,seeds); % set locations of seeds to be -Inf (cause matlab watershed)
%figure('name','cyto_min','NumberTitle', 'off')
%imshow(cyto_min,[]);
cyto_ws = watershed(cyto_min);
cyto_ws(cyto_mask==0)=0; % remove areas that aren't in our cyto mask
figure('name','cyto_ws','NumberTitle', 'off')
imshow(cyto_ws,[]);

labelled_cyto = cyto_ws; % set a better variable name

% Clear cells touching the boarder
labelled_cyto = imclearborder(labelled_cyto);

% Fill holes
labelled_cyto = imfill(labelled_cyto,'holes');
figure('name','labelled_cyto','NumberTitle', 'off')
imshow(labelled_cyto,[]);



%% Color overlay v1
% Display original image
figure
imshow(imoverlay(uint8(cyto*255),bwperim(labelled_cyto),'r'),[]);
hold on
% Display color overlay
labelled_cyto_rgb = label2rgb(uint32(labelled_cyto), 'jet', [1 1 1], 'shuffle');
himage = imshow(labelled_cyto_rgb,[]);
himage.AlphaData = cyto_mask*.3;
% Display red borders
[xm,ym]=find(seeds);
hold on
plot(ym,xm,'or','markersize',2,'markerfacecolor','g','markeredgecolor','g')

%% Color overlay v2 (perim only)
% Display original image
figure
imshow(uint8(cyto*255),[]);
hold on
% Display color overlay
labelled_cyto_rgb = label2rgb(uint32(bwlabel(bwperim(labelled_cyto))), 'jet', [1 1 1], 'shuffle');
himage = imshow(labelled_cyto_rgb,[]);
himage.AlphaData = cyto_mask*.3;
% Display red borders
[xm,ym]=find(seeds);
hold on
plot(ym,xm,'or','markersize',2,'markerfacecolor','g','markeredgecolor','g')