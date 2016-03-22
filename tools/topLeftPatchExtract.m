function D = topLeftPatchExtract(im,patchSiz)

[m,n,~] = size(im);
C = padarray(im,[patchSiz-1,patchSiz-1],'replicate','post');
A = im2col(C, [patchSiz,patchSiz], 'sliding');
A = A';
D = reshape(A,[m,n,patchSiz*patchSiz]);
D = uint8(D*255);
