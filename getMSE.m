function MSE = getMSE(img1,img2)

if any(size(img1) ~= size(img2))
    nRows = min(size(img1,1),size(img2,1));
    nCols = min(size(img1,2),size(img2,2));
    img1 = img1(1:nRows,1:nCols,:);
    img2 = img2(1:nRows,1:nCols,:);
    %return;
end

% calcule o erro
err = double(img1)-double(img2);

% calcule o MSE
MSE = sum(err(:).^2)/numel(err);

end