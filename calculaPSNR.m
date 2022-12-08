function Valor_psnr = calculaPSNR(A,B)
% PSNR (Peak Signal to noise ratio)

if (size(A) ~= size(B))
   error('O tamanho das matrizes s�o diferentes!')

   Valor_psnr = NaN;
   return; 
elseif (A == B)
   disp('As imagens s�o identicas: a PSNR tem valor infinito')

   Valor_psnr = Inf;
   return;   
else

    maxValue = double(max(A(:)));

    % C�lculo MSE.
    mseImage = (double(A) - double(B)) .^ 2;
    [rows, columns] = size(A);
    
    mse = sum(mseImage(:)) / (rows * columns);

    % C�lculo da PSNR 
    Valor_psnr = 10 * log10( 256^2 / mse);
end


end

