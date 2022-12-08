function varargout = orthogen2(alpha,m);

% h = orthogen2(alpha,m)
% [h,g] = orthogen2(alpha,m)
% [h,g,alpha_new] = orthogen2(alpha,m)
%
% Constructs an array, h(1:2N) of lowpass orthonormal
% FIR filter coefficients for any N. The filter can have
% 0, 1 or 2 vanishing moments.
%
% g(1:2N) is an array of highpass filter coefficients, obtained by 
% alternating flip of h. 
%
% Inputs:
% - alpha: array (1...N) with N free parameters which are angles in
%          radians
% - m: number of desired vanishing moments (0, 1 or 2)
%      (default: m=0)
%
% If one or two vanishing moments are required, alpha(N) is replaced
% by a value which ensure the first vanishing moment.
%
% If two vanishing moments are required, alpha(N-1) is also replaced by a
% convenient value. The solution for alpha(N-1) is not guaranteed to exist
% and depends on the choice of alpha(1:N-2). If the solution does not exist
% for the provided set of angles, then alpha(N-2), alpha(N-3),..., are also
% replaced until the existence of the solution is achieved.
%
% Output alpha_new returns the value of alpha with the changed values.
%
% This function is based on function orthogen.m, presented in 
% Sherlock & Monro, IEEE Trans. Signal Proc., 46 (1998) 1716-1720.
%
% Henrique M. Paiva, Oct 2007

% input verification
  if ~exist('m','var'), m = 0; end

% adjust alpha in order to ensure the vanishing moments
  switch m  
      case 1, 
          alpha_new = firstMoment(alpha);
      case 2,
          alpha_new = secondMoment(alpha);
          alpha_new = firstMoment(alpha_new);
      otherwise,
          alpha_new = alpha;
  end

% call orthogen
  h = orthogen(alpha_new);

% obtain g by alternating flip of h  
  g = get_g(h);
  
% set outputs 
  if nargout > 0, varargout{1} = h; end
  if nargout > 1, varargout{2} = g; end
  if nargout > 2, varargout{3} = alpha_new; end

return

%-------------------------------------------------------------
% Subfunction get_g
% - Obtain g by alternating flip of h
%-------------------------------------------------------------
function g = get_g(h)
g = h(end:-1:1);
i = 1:2:length(g);
g(i) = -g(i);
return

%-------------------------------------------------------------
% Subfunction firstMoment
% - Change alpha in order to ensure one vanishing moment
%-------------------------------------------------------------
function alpha_new = firstMoment(alpha)
alpha_new = alpha;
alpha_new(end) = pi/4 - sum(alpha_new(1:end-1));
return

%-------------------------------------------------------------
% Subfunction secondMoment
% - Change alpha in order to ensure two vanishing moments
% - Note: Function firstMoment shall be called after secondMoment
%  . e.g: alpha_new = secondMoment(alpha); alpha_new = firstMoment(alpha_new);
% - Required internal functions: get_sum and get_mu
%-------------------------------------------------------------
function alpha_new = secondMoment(alpha)

% Get number of elements 
N = length(alpha);

% Check if solution in alpha(N-1) exists. If not, change alpha(N-2),
% alpha(N-3),..., until the existence of the solution is achieved 

% lower bound check
j = N-2;
while get_sum(alpha) < -3/2 
    for el = j:N-2
        alpha(el) = 0.5*(pi/2-2*sum(alpha(1:el-1)));
    end    
    j = j-1;
end    
% higher bound check
j = N-2;
while get_sum(alpha) > 1/2 
    for el = j:N-2
        alpha(el) = 0.5*(-pi/2-2*sum(alpha(1:el-1)));
    end    
    j = j-1; 
end    

% replace alpha(N-1) by the value which ensures the second vanishing moment
alpha(N-1) = (1/2)*asin(-1/2-get_sum(alpha))-sum(alpha(1:N-2));

% set output
alpha_new = alpha;

return


%-------------------------------------------------------------
% Subfunction get_sum
% - Return the following sum:
%   sum{k=1 to N-2} [sin(sum{i=1 to k}{2*alpha(i)} ]
% - Note: mu(k) = sum{i=1 to k}{2*alpha(i)}
%-------------------------------------------------------------
function sum1 = get_sum(alpha)
N = length(alpha);
mu = get_mu(alpha);
sum1 = sum(sin(mu(1:N-2)));
return

%-------------------------------------------------------------
% Subfunction get_mu
% - Return the following sum:
%   mu(k) = sum{i=1 to k}{2*alpha(i)}
%-------------------------------------------------------------
function mu = get_mu(alpha)
N = length(alpha);
mu = zeros(1,N);
for k = 1:N
    mu(k) = 2*sum(alpha(1:k));
end    
return

%-------------------------------------------------------------
% Subfunction orthogen
%-------------------------------------------------------------

function h = orthogen(alpha);
%h = orthogen(alpha)
% Constructs an array, h(1...N) of lowpass orthonormal
% FIR filter coefficients for any even N >= 2.
%
% The input array, alpha(1...N/2) gives N/2 free
% parameters which are angles in radians.  If the
% angles sum to pi/4 the filter corresponds to a 
% regular wavelet.
%
% Sherlock & Monro, IEEE Trans. Signal Proc., 46 (1998) 1716-1720.


N  = 2*length(alpha);
h  = zeros(1, N);
lo = N/2;
hi = lo + 1;
h(lo) = cos(alpha(1));
h(hi) = sin(alpha(1));
nstages = N/2;


for stage = 1 : nstages-1

	c = cos(alpha(stage+1));
	s = sin(alpha(stage+1));
		
	h(lo-1) = c*h(lo);
	h(lo)   = s*h(lo);
	h(hi+1) = c*h(hi);
	h(hi)   = -s*h(hi);
	
	nbutterflies = stage-1;
	butterflybase = lo+1;

	for butterfly = 1 : nbutterflies
		hlo = h(butterflybase);	
		hhi = h(butterflybase+1);
		h(butterflybase)   = c*hhi - s*hlo;
		h(butterflybase+1) = s*hhi + c*hlo;
		butterflybase = butterflybase + 2;
    end
	
	lo = lo - 1;
	hi = hi + 1;
    
end
	