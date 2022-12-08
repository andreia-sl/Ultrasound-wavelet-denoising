function compareFilters

% lista de wavelets para comparacao
wavList = {
    {'bior3.3', 'db4', 'sym4'} % 8 taps
    {'bior1.5', 'bior2.4','bior4.4', 'db5', 'sym5'} % 10 taps
    {'bior3.5', 'bior5.5', 'db6', 'sym6', 'coif2'} % 12 taps
    {'bior3.7', 'db8', 'sym8'} % 16 taps
    {'db10','sym10'}; % 20 taps
    };

% lista de diretorios
inputDirList = {
    'Appendix'
    'Gallbladder'
    'Gastrointestinal tract'
    'Liver'
    'Pancreas'
    'Spleen'
    };

% lista de arquivos
fileList = cell(length(inputDirList),1);
for iDir = 1:length(inputDirList)
    currDir = fullfile(cd,'UltraSoundCases\',inputDirList{iDir});
    currFileList = dir(fullfile(currDir,'*.jpg'));
    currFileList = arrayfun(@(f) fullfile(currDir,f.name),currFileList, 'UniformOutput', false);
    fileList{iDir} = currFileList;
end

% lista de arquivos para validacao
fileList_val = cellfun(@(f) f(2:end)',fileList,'UniformOutput',false);
fileList_val = [fileList_val{:}]';

% flags
flags = struct;
flags.HardThresh        = 0;  % flag: 0 - soft thresholding; 1 - hard thresholding
flags.OptThresh         = 0;  % otimizar tambem o limiar
flags.OneFilterPerLevel = 0;  % flag: um filtro diferente para cada nivel
flags.defaultWavTree    = 1;  % flag: usar decomposicao por wavelet tree padrao

noiseVarArray = 0.02:0.02:0.2;
if 1
    noiseVarArray = linspace(0.01,0.2,7); warning('Reduzindo os valores de ruido para teste');
end

% dwt mode
dwtmode('per');
nLevels = 3;

% para cada numero de taps
for idxTaps = 1:size(wavList,1)
    
    % wavelets para comparar
    waveletsToCompare = wavList{idxTaps};
    
    % quantidades
    nFile  = length(fileList_val);
    nNoise = length(noiseVarArray);
    nWav   = length(waveletsToCompare);
    % inicializacoes
    PSNR_array_NoFilt = zeros(nFile,nNoise);
    PSNR_array_Tta0   = zeros(nFile,nNoise,nWav);
    SSIM_array_NoFilt = zeros(nFile,nNoise);
    SSIM_array_Tta0   = zeros(nFile,nNoise,nWav);
    RMSE_array_NoFilt = zeros(nFile,nNoise);
    RMSE_array_Tta0   = zeros(nFile,nNoise,nWav);
    
    % numero de casos
    nCases = nFile*nNoise;
    loopStart = now;
    hBar = waitbar(0);
    
    % loop de analise
    for iFile = 1:nFile
        currFile = fileList_val{iFile};
        currImg = imread(currFile);
        %keyboard
        for iNoise = 1:nNoise
            % informacoes sobre o loop
            iCase = (iFile-1)*nNoise + iNoise;
            loopEnd = loopStart + (nCases/iCase)*(now-loopStart);
            infoStr = sprintf('Case %i of %i\nStart: %s, End: %s',iCase,nCases,...
                datestr(loopStart,'HH:MM:SS'),datestr(loopEnd,'HH:MM:SS'));
            waitbar(iCase/nCases,hBar,infoStr);
            % imagem atual com ruido
            currImg = imresize(currImg,[200 200]);
            currImg_Noise = imnoise(currImg,'speckle',noiseVarArray(iNoise));
            % calculo sem processamento
            currPSNR_noFilt = calculaPSNR(currImg,currImg_Noise);
            currSSIM_noFilt = ssim(currImg,currImg_Noise);
            currRMSE_noFilt = sqrt(getMSE(currImg,currImg_Noise));
            % calculo com wavelets padrao
            for iWav = 1:nWav
                currWavelet = waveletsToCompare{iWav};
                h = wfilters(currWavelet);
                currTheta = parameterize2(h);
                currTheta = currTheta(1:end-1);
                if flags.OneFilterPerLevel
                    currTheta = repmat(currTheta(:),nLevels,1);
                end
                [~,currPSNR_Tta0,~,currSSIM_Tta0,currRMSE_Tta0]   = mycost(currTheta, currImg,currImg_Noise,nLevels,flags);
                PSNR_array_Tta0(iFile,iNoise,iWav)   = currPSNR_Tta0;
                SSIM_array_Tta0(iFile,iNoise,iWav)   = currSSIM_Tta0;
                RMSE_array_Tta0(iFile,iNoise,iWav)    = currRMSE_Tta0;
            end
            % armazene os resultados
            PSNR_array_NoFilt(iFile,iNoise) = currPSNR_noFilt;
            %
            SSIM_array_NoFilt(iFile,iNoise) = currSSIM_noFilt;
            %
            RMSE_array_NoFilt(iFile,iNoise) = currRMSE_noFilt;
        end
    end
    delete(hBar);
    
    legStr = {waveletsToCompare{:},'No Filter'};
    
    figure
    set(gcf,'units','pixels','pos',[100 60 550 930]);
    
    % PLOT DOS VALORES MEDIOS DE PSNR
    subplot(3,1,1);
    hold all; grid on; box on;
    for iWav = 1:nWav
        tHdl = smoothplot(noiseVarArray,mean(PSNR_array_Tta0(:,:,iWav),1),'-.','linewidth',2);
    end
    smoothplot(noiseVarArray,mean(PSNR_array_NoFilt,1),'--','linewidth',2);
    legend(legStr);
    xlabel('Noise Variance');
    ylabel('PSNR');
    
    
    % PLOT DOS VALORES MEDIOS DE SSIM
    subplot(3,1,2);
    hold all; grid on; box on;
    for iWav = 1:nWav
        tHdl = smoothplot(noiseVarArray,mean(SSIM_array_Tta0(:,:,iWav),1),'-.','linewidth',2);
    end
    smoothplot(noiseVarArray,mean(SSIM_array_NoFilt,1),'--','linewidth',2);
    legend(legStr);
    xlabel('Noise Variance');
    ylabel('SSIM');
    
    % PLOT DOS VALORES MEDIOS DE RMSE
    subplot(3,1,3);
    hold all; grid on; box on;
    for iWav = 1:nWav
        tHdl = smoothplot(noiseVarArray,mean(RMSE_array_Tta0(:,:,iWav),1),'-.','linewidth',2);
    end
    smoothplot(noiseVarArray,mean(RMSE_array_NoFilt,1),'--','linewidth',2);
    legend(legStr);
    xlabel('Noise Variance');
    ylabel('RMSE');

end

keyboard

return

function tHdl = smoothplot(x,y,varargin)
    xr = linspace(min(x),max(x),100);
    yr = interp1(x,y,xr,'spline');
    tHdl = plot(xr,yr,varargin{:});
return