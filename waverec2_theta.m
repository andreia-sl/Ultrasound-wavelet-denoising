function wavCoefRec = waverec2_theta(wavCoef,theta,flagDefaultWavTree)

    % n�mero de n�veis de decomposi��o
    n = size(wavCoef,1)-1;

    % Determina o n�mero dos N �ngulos para cada filtro passa-baixa
    N = length(theta)/n;

    % Inicializa��o da sa�da
    wavCoefRec = cell(size(wavCoef));

    % Copie todos os coeficientes do ultimo nivel
    wavCoefRec(end,:) = wavCoef(end,:);
    % se arvore wavelet padrao, copie tambem os coeficientes dos canais
    % passa-altas
    if flagDefaultWavTree
        wavCoefRec(2:end,2:4) = wavCoef(2:end,2:4);
    end

    % for each decomposition level
    for jLevel = n:-1:1

        % Step 1: obter os N �ngulos do filtro para o n�vel atual
        thetak = theta((1 + (jLevel-1)*N) : jLevel*N);

        % Step 2: obter os filtros wavelet utilizando
        thetak(end+1) = Inf; % inicializa com 'Inf' (ser� sobrescrito em orthogen2)
        [h,g] = orthogen2(thetak,1);
        hr = h(end:-1:1);
        gr = g(end:-1:1);

        % for each parent node
        if flagDefaultWavTree
            nodesToReconstruct = 1; % decomponha apenas o canal passa-baixas
        else
            nodesToReconstruct = 1:4^(jLevel-1); % decomponha todos os canais
        end
        for iNode = nodesToReconstruct
            % leia os coeficientes
            CA = wavCoefRec{jLevel+1,4*(iNode-1)+1};
            CH = wavCoefRec{jLevel+1,4*(iNode-1)+2};
            CV = wavCoefRec{jLevel+1,4*(iNode-1)+3};
            CD = wavCoefRec{jLevel+1,4*(iNode-1)+4};
            % adjust dimensions
            if numel(CA) ~= numel(CH)
                CA = CA(1:size(CH,1),1:size(CH,2),1:size(CH,3));
            end
            % processando
            currOutput = idwt2(CA,CH,CV,CD,hr,gr);
            % currOutput
            wavCoefRec{jLevel,iNode} = currOutput;
        end
        
end

