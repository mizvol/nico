function [Ahat,pihat] = nico(X,n,APri,piPri)
% NICO   Network Inference from Co-Occurrences
%
%   Calling syntax:
%      [Ahat,pihat] = nico(X,N,APri,piPri)
%
%   Inputs:
%      X - Cell array of co-occurrence data.  Each cell contains one
%          unordered list (column vector) of nodes in a path
%      N - Number of nodes in the network
%      APri - Prior on the transition matrix (Dirichlet parameters)
%      piPri - Prior on the initial state distribution (Dirichlet params)
%   Output:
%      Ahat - Estimated transition matrix
%      pihat - Estimated initial state distribution

% Michael Rabbat
% McGill University
% michael.rabbat@mcgill.ca
%
% 15 June 2007 - New version in Matlab code using sparse matrices

%% Argument checking
if (nargin < 4)
    error(['Not enough input arguments.  ' ...
        'If you do not have priors, enter zeros']);
end

%MGR Add more error checking here later

%% Basic info about the data
% Number of paths
T = size(X,1);
% Number of nodes in each path
Nm = zeros(T,1);
for m=1:T
    Nm(m) = length(X{m});
end

%% Initialize pihat and Ahat
% Assume all states appear at least once in the data
pihat = 1 + 0.3*rand(n,1);
pihat = pihat./sum(pihat);

% Construct Ahat as a sparse matrix
% First determine an upperbound on the number of non-zero entries
maxnnz = 0;
for m=1:T
	maxnnz = maxnnz + nchoosek(Nm(m),2);
end

% Determine initial sparse matrix structure based on the data
ii = zeros(maxnnz,1);
jj = zeros(maxnnz,1);
ss = ones(maxnnz,1);
last = 0;
for m=1:T
    V = nchoosek(X{m},2);
	numadd = size(V,1);
    ii(last+1:last+numadd) = V(:,1);
	jj(last+1:last+numadd) = V(:,2);
	last = last + numadd;
end
Ahat = sparse(ii,jj,ss,n,n);
Ahat = (Ahat + Ahat')./2;
Ahat = spones(Ahat) + 0.4*sprand(Ahat);
clear ii jj ss;

% Normalize rows
for i=1:n
	rowsum = sum(Ahat(i,:));
	if (rowsum > 0)
		Ahat(i,:) = Ahat(i,:)./rowsum;
	end
end

% Note to self: For normalizing over rows of a sparse matrix it's probably
% more efficient to store the transpose of the adjacency matrix...

%% EM Algorithm Iterations
tol = 0.01;
Kmax = 100;
plcutoff = 10;
RAlpha = cell(T,2);
for k=1:Kmax
	
	% E-Step
	for m=1:T
		gamma = zeros(Nm(m),1);
		Gamma = zeros(Nm(m),Nm(m));
		if ((Nm(m) > 1) && (Nm(m) <= plcutoff))
			% Exact E-Step
			tau = zeros(Nm(m),1);
			gamma = zeros(Nm(m),1);
			Gamma = zeros(Nm(m),Nm(m));
			for t=1:Nm(m)
				[gamma,Gamma] = permrecurse(t,1,Nm(m),tau,X{m},gamma,...
					Gamma,pihat,Ahat);
			end
			
		elseif (Nm(m) > plcutoff)
			
			% Approximate E-Step
			numSamples = max(10000, 2*(Nm(m))^4);
			
			for i=1:numSamples
				%% Sample a permutation
				tau = zeros(Nm(m),1);
				% Initialize flag variable
				f = ones(Nm(m),1);
				% Sample the first slot in the path
				pprime = pihat(X{m});
				cumprobs = cumsum(pprime./sum(pprime));
				larger = find(cumprobs >= rand);
				s = larger(1);
				f(s) = 0;
				tau(1) = s;
				w = 1;
				Aprime = Ahat(X{m},X{m});
				% Sample the remaining slots
				for t=2:Nm(m)
					pprime = Aprime(s,:)'.*f;
					if (sum(pprime) == 0)
						break;
					end
					w = w*sum(pprime);
					cumprobs = cumsum(pprime./sum(pprime));
					larger = find(cumprobs >= rand);
					s = larger(1);
					tau(t) = s;
					f(s) = 0;
				end
				
				% Bookkeeping
				if (length(find(tau==0)) == 0)
					gamma(tau(1)) = gamma(tau(1)) + w;
					for t=2:Nm(m)
						Gamma(tau(t-1),tau(t)) = Gamma(tau(t-1),tau(t)) + w;
					end
				else
					% Didn't sample a complete permutation
				end
			end
			
		else
			% Path only has one node
		end
		
		% Store info for this path
		RAlpha{m,1} = gamma./sum(gamma);
		RAlpha{m,2} = Gamma./sum(gamma);
	end
	
	% M-Step
	C = APri;
	c = piPri;
	for m=1:T
		for t=1:Nm(m)
			c(X{m}(t)) = c(X{m}(t)) + RAlpha{m,1}(t);
		end
		for t=1:Nm(m)-1
			for tt=t+1:Nm(m)
				C(X{m}(t),X{m}(tt)) = C(X{m}(t),X{m}(tt)) + RAlpha{m,2}(t,tt);
				C(X{m}(tt),X{m}(t)) = C(X{m}(tt),X{m}(t)) + RAlpha{m,2}(tt,t);
			end
		end
	end
	AhatOld = Ahat;
	pihatOld = pihat;
	Ahat = max(C,0);
	pihat = max(c,0);

	% Normalize
	if (sum(pihat)<=0)
		% Hmm, something went wrong
		error('All entries of pihat vanished');
	end
	pihat = pihat./sum(pihat);
	for i=1:n
		s = sum(Ahat(i,:));
		if (s > 0)
			Ahat(i,:) = Ahat(i,:)./s;
		end
	end
	
	% Compute change in Q
	Q = 0;
	Qold = 0;
	for m=1:T
		for t=1:Nm(m)
			Q = Q + RAlpha{m,1}(t)*log(pihat(X{m}(t)));
			Qold = Qold + RAlpha{m,1}(t)*log(pihatOld(X{m}(t)));
		end
		for t=1:Nm(m)-1
			for tt=t+1:Nm(m)
				Q = Q + RAlpha{m,2}(t,tt)*log(Ahat(X{m}(t),X{m}(tt)) + eps);
				Q = Q + RAlpha{m,2}(tt,t)*log(Ahat(X{m}(tt),X{m}(t)) + eps);

				Qold = Qold + RAlpha{m,2}(t,tt)*log(AhatOld(X{m}(t),X{m}(tt)) + eps);
				Qold = Qold + RAlpha{m,2}(tt,t)*log(AhatOld(X{m}(tt),X{m}(t)) + eps);
			end
		end
	end
	Delta = Q - Qold;
	
	% Check stopping criterion
	%sc = norm(AhatOld - Ahat, 'fro')/tol;
	sc = Delta/tol;
	%disp(num2str(sc));
	fprintf(1,'Iter=%4d, Delta=%5.4f, Q=%5.4f\n', ...
		k, Delta, Q);
	if (sc < 1)
		break;
	end
end

if ((k == Kmax) && (sc > 1))
	disp(['Warning: Number of EM iterations exceeded max.  Stopped after '...
		num2str(k) ' iterations']);
else
	disp(['Terminated successfully after ' num2str(k) ' iterations']);
end

return;


%% Recursive evaluation of all permutations in exact E-step
function [gamma, Gamma] = permrecurse(visit, level, Nm, tau, x, gamma, ...
	Gamma, pihat, Ahat)

tau(visit) = level;
if (level == Nm)
	% Compute likelihood of this link ordering
	p = pihat(x(tau(1)));
	for t=2:Nm
		p = p*Ahat(x(tau(t-1)),x(tau(t)));
	end
	
	% Bookkeeping
	gamma(tau(1)) = gamma(tau(1)) + p;
	for t=2:Nm
		Gamma(tau(t-1),tau(t)) = Gamma(tau(t-1),tau(t)) + p;
	end
else
	% Continue recursive DF traversal down the tree
	for t=1:Nm
		if (tau(t)==0)
			[gamma, Gamma] = permrecurse(t, level+1, Nm, tau, x, gamma, ...
				Gamma, pihat, Ahat);
		end
	end
end
tau(visit) = 0;

return;
