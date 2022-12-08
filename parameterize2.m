function alpha = parameterize2(h);

% alpha = parameterize2(h)
% h(1...N) - array of lowpass orthonormal FIR filter coefficients, N even >= 2.
% The output array alpha(1...N/2) consists of N/2 angles in radians
% that parameterize filter h using the formulation presented by
% Sherlock & Monro, IEEE Trans. Signal Proc., 46 (1998) 1716-1720.
% This function is intended to be the inverse of function orthogen.m, 
% in the sense that, if 
%     alpha2 = parameterize2(orthogen(alpha))
% then
%     orthogen(alpha2) = orthogen(alpha)          
% Note that alpha2 is not necessarily equal to alpha, 
% because function orthogen.m is not injective. 
% Author: Henrique M. Paiva, April 30th, 2003
% Latest version: July 18th, 2003

% Number of elements of h
  N = length(h);

% Input consistency verification
  if((N/2)~=floor(N/2))
      Error('h must have an even number of elements');
  end

% Indices over array h
  lo = 1;
  hi = N;

% Number of stages
  nstages = N/2;

% Initialization   
  alpha = zeros(1,nstages);

% Determination of the angular parameters alpha

  for stage = nstages:-1:1

    % Determination of current angle alpha  
      alpha(stage) = atan2(h(lo+1),h(lo));

    % Determination of filter (k-1) as a function of filter (k)      
      % flow diagram weights
        c = cos(alpha(stage));
        s = sin(alpha(stage));
      % external weights
        if (s~=0)
          h(lo+1) =  h(lo+1)/s;
          h(hi-1) = -h(hi-1)/s;
        else  
          h(lo+1) =  h(lo)/c;
          h(hi-1) =  h(hi)/c;
        end
      % internal parts (butterflies)
  	    nbutterflies  = stage-2;
   	    butterflybase = lo+2;
	    for butterfly = 1 : nbutterflies
  		   hlo = h(butterflybase);	
		   hhi = h(butterflybase+1);
		   h(butterflybase)   = c*hhi - s*hlo; %internal weight
		   h(butterflybase+1) = s*hhi + c*hlo; %internal weight
		   butterflybase = butterflybase + 2;
           % Note: The inverse butterfly is identical to the direct butterfly.
        end
    % Index change    
   	  lo = lo + 1;
	  hi = hi - 1;
  end

