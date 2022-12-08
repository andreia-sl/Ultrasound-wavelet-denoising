function [J,PSNR,NewImgh,SSIM,RMSE] = mycost(theta,img,imgR,nLevels,flags)

    % decomposicao das entradas
    if flags.OptThresh
        thresh_fac = theta(end);
        theta = theta(1:end-1);
    else
        thresh_fac = 1; % mantenha o limiar padrao
    end

    % repetir theta (usar o mesmo em todos os niveis)
    if ~flags.OneFilterPerLevel
        theta = repmat(theta(:),nLevels,1);
    end

    % decomposicao wavelet
    % dwtmode('per');
    wavCoef = wavedec2_theta(imgR,nLevels,theta,flags.defaultWavTree);

    % defina os indices dos nos na arvore
    if ~flags.defaultWavTree
        nodeIdx = [repmat(nLevels,4^nLevels,1) (1:4^nLevels)']; % arvore wavelet completa
    else
        [a,b] = meshgrid(2:nLevels,2:4); 
        nodeIdx = [a(:) b(:)];
        nodeIdx(end+1,1:2) = [nLevels 1];
    end

    % coeficientes dos nos folha
    nodeCoef = cell(size(nodeIdx,1),1);
    for i = 1:size(nodeIdx,1)
        nodeCoef{i} = (wavCoef{nodeIdx(i,1),nodeIdx(i,2)}(:))';
    end
    nodeCoef = [nodeCoef{:}];

    % valor default de threshold
    sigma = median(abs(nodeCoef))/0.6745;
    thresh = sigma*sqrt(2*log(numel(imgR))); % T é o valor do limiar (threshold)

    % multiplique pela entrada thresh_fac
    thresh = thresh*thresh_fac;

    % Thresholding - Soft ou Hard Thresholding
    % (aplicado aos coeficientes do ultimo nivel)
    for i = 1:size(nodeIdx,1)
        i1 = nodeIdx(i,1);
        i2 = nodeIdx(i,2);
        wavCoef{i1,i2}(abs(wavCoef{i1,i2})<thresh) = 0;
        if ~flags.HardThresh
            wavCoef{i1,i2}(wavCoef{i1,i2}<-thresh) = wavCoef{i1,i2}(wavCoef{i1,i2}<-thresh)+thresh;
            wavCoef{i1,i2}(wavCoef{i1,i2}> thresh) = wavCoef{i1,i2}(wavCoef{i1,i2}>thresh) -thresh;
        end
    end

    % recomposicao
    wavCoefRec = waverec2_theta(wavCoef,theta,flags.defaultWavTree);

    % construa NewImgh 
    NewImgh = wavCoefRec{1,1};

    % transforme a imagem em tipo inteiro
    NewImgh = cast(NewImgh,'like',img);

    % Testar com duas psnr's de imagens contaminadas com ruídos diferentes
    PSNR=calculaPSNR(NewImgh,img);

    SSIM = ssim(NewImgh,img);

    MSE = getMSE(NewImgh,img);
    RMSE = sqrt(MSE);

J = -PSNR;
