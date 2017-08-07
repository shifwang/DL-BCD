% Generate a toy dictionary learning problem

% underlying dictionary: dict
clc, clear all, close all;
width         = 32;
height        = 32;
dict_size     = 9;
square_length = 8;
square_width  = 3;
add_noise     = false; % true if add Gaussian noise with deviation 0.1
dict          = zeros(width*height, dict_size);
for i = 1:dict_size
    angle = pi / dict_size * i;
    img   = zeros(width, height);
    for x = 1:width
        for y = 1:height
            tmp1 = cos(angle)*(x - width/2) + sin(angle)*(y - height/2);
            tmp2 = -sin(angle)*(x - width/2) + cos(angle)*(y - height/2);
            if abs(tmp1) < square_length && abs(tmp2) < square_width
                img(x,y) = 1;
            end
        end
    end
    dict(:,i) = reshape(img,width*height,1);
end

%% Visualize the dictionary
figure;
for i = 1:dict_size
    subplot(3, 3, i)
    imshow(reshape(dict(:,i),width, height));
end

%% Generate samples
positive       = false;
sparsity       = 5;
distribution   = 'Gaussian';
sample_size    = 1000;
if strcmp(distribution,'Gaussian')
    coef = zeros(dict_size,sample_size);
    if positive
        for col = 1:sample_size
            coef(randperm(dict_size,sparsity),col) = abs(randn(sparsity,1));
        end
    else
        for col = 1:sample_size
            coef(randperm(dict_size,sparsity),col) = randn(sparsity,1);
        end
    end
end
if add_noise
    data = dict * coef + randn(width * height, sample_size)*0.1;
else
    data = dict * coef;
end

%% Visualize some samples
figure;
shading interp;
axis off;
for i = 1:12
    min_value = min(data(:,i));
    max_value = max(data(:,i));
    subplot(4, 3, i)
    imshow(reshape((data(:,i)-min_value)./(max_value - min_value),width, height));
end

