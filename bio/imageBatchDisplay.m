function imageBatchDisplay(D,width,height,nrow,ncol,colorScale,r)
% D: a list of images -- each column corresponding to an image
% width/height: the width/height of an image
% nrow/ncol: the number of rows/cols in the layout?
% colorScale: can be '' or 'gray'
% r: (optional) a scaling parameter
   [p,n] = size(D);
   nR = nrow*height;
   nC = ncol*width;
   I = zeros(nR,nC);
   rowStart = 1;
   colStart = 1;
   k = 1;
   for i = 1:nrow
       for j = 1:ncol
           indR = (rowStart + (i-1)*height):( i*height);
           indC = (colStart + (j-1)*width):( j*width);
           I(indR,indC) = reshape(D(:,k),height,width);
           k = k+1;
           if k>n
               break;
           end
       end
       if k>n
           break;
       end
   end
   %if ~exist('r')
   %    r = [min(min(I)),max(max(I))];
   %end
   figure('position',[1404,0,150*ncol,75*nrow]);
   if strmatch(colorScale,'gray','exact')
       if ~exist('r')
           r = [min(min(1-I)),max(max(1-I))];
       end
       if min(r) == max(r)
           imagesc(1-I);
       else
           imagesc(1-I,r);
       end
       colormap(gray);
   else
       if ~exist('r')
           if false
               r = [min(min(I)),quantile(I(:),.999)];
           else
               top = quantile(abs(I(:)),.999);
               r = [-top, top];
       end
       if min(r) == max(r)
           imagesc(I);
       else
           %I = (I - min(r))/max(r);
           colormap('redblue')
           imagesc(-I, r);
       end

   end
   
   axis off;        
  
   for j = 1:ncol
       x = [colStart + (j-1)*width,1];
       y = [colStart + (j-1)*width,nR];
       drawLineSeg(x,y);
   end
   for i = 1:nrow
       x = [1,rowStart + (i-1)*height];
       y = [nC,rowStart + (i-1)*height];
       drawLineSeg(x,y);
   end
   
   k = 1;
   for i = 1:nrow
       for j = 1:ncol
           y =  i*height - 8;
           x =  j*width - 8;
           text(x,y,num2str(k),'color','red','horizontalAlignment','center');
           k = k+1;  
           if k > n
               break;
           end
       end
       if k > n
           break;
       end
   end
end
