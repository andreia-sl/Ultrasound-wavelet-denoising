function wavCoef = wavedec2_theta(x,n,theta,flagDefaultWavTree)

    % Determine o número N de ângulos para cada filtro passa-baixa
    N = length(theta)/n;

    % Inicialização da saída
    if flagDefaultWavTree
        wavCoef = cell(n+1,4);
    else
        wavCoef = cell(n+1,4^n);
    end
    wavCoef{1,1} = x; % wavCoef{1,1} contém os dados de entrada

    % para cada nível de decomposição
    for jLevel = 1:n

        % Step 1: obter os N ângulos do filtro para o nível atual
        thetak = theta(1 + (jLevel-1)*N : jLevel*N);

        % Step 2: obter os filtros wavelet utilizando orthogen
        thetak(end+1) = Inf; % inicialuza cin 'Inf' (será sobescrito em orthogen2)
        [h,g] = orthogen2(thetak,1);

        % para cada nó pai
        if flagDefaultWavTree
            nodesToDecompose = 1; % decomponha apenas o canal passa-baixas
        else
            nodesToDecompose = 1:4^(jLevel-1); % decomponha todos os canais
        end
        for iNode = nodesToDecompose
            % curr input
            currInput = wavCoef{jLevel,iNode};
            % processando
            [CA,CH,CV,CD] = dwt2(currInput,h,g);
            % salva as saídas
            wavCoef{jLevel+1,4*(iNode-1)+1} = CA;
            wavCoef{jLevel+1,4*(iNode-1)+2} = CH;
            wavCoef{jLevel+1,4*(iNode-1)+3} = CV;
            wavCoef{jLevel+1,4*(iNode-1)+4} = CD;
        end
        
end

