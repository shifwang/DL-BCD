function imageBatchDisplaySaveMemory(D,width,height,nrow,ncol,colorScale,r, template)
% D: a list of images -- each column corresponding to an image
% width/height: the width/height of an image
% nrow/ncol: how many rows/cols in the layout?
% colorScale: can be '' or 'gray'
% r: a scaling parameter, 0.5 for 64by32 image and 0.25 for 32by16 images
% e.g. imageBatchDisplaySaveMemory(D,32,16,8,8,'gray',0.25);
% template = generateTemplate();

%if ~exist('r')
%    template = imresize(template,0.5,'nearest');
%else
%    template = imresize(template,r,'nearest');
%end
%for i = 1:size(D,2)
%    temp = D(:,i);
%    maxAbs = max(abs(temp));
%    if max(temp) == maxAbs
%       D(:,i) = temp/max(temp);
%    else
%       D(:,i) = -temp/maxAbs;
%    end
%end

ind = find(template(:,:,1)==1);
Dfull = zeros(height*width,size(D,2));
Dfull(ind,:) = D;
imageBatchDisplay(Dfull, width, height, nrow,ncol,colorScale);
%add boundary
boundary = [];
for i = 1:size(template, 1)
    for j = 1:size(template, 2)
        if template(i, j, 1) == 1 && (i == 1 || template(i - 1, j, 1) == 0 || i == size(template, 1) || template(i + 1, j, 1) == 0 || j == 1 || j == size(template, 2))
            boundary = [boundary; [i, j]];
        end
    end
end
hold on;
thres = 1;
for i = 1:nrow
    for j = 1:ncol
        if thres > size(D, 2)
            break
        end
        h = scatter((j - 1) * width + boundary(:, 2), (i - 1) * height + boundary(:, 1), 'filled','k', 'SizeData', 20 );
        alpha(h, .3)
        thres = thres + 1;
    end
    if thres > size(D, 2)
        break
    end
end