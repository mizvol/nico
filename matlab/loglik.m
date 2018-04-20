function ll = loglik(X,A,Pi);
% LOGLIK   Compute the NICO log likelihood
%
%    loglik(X,A,Pi) computes P(X | A, Pi) for the NICO model.
%    For paths with no more than 10 nodes this function explicitly computes
%    the log-likelihood term (enumerating all permutations).  For longer
%    paths a Monte Carlo estimate is used.
%
%    Inputs:
%       X - a cell array of shuffled observations
%       A - estimated transition matrix
%       Pi - estimated initial state distribution

% Michael Rabbat
% McGill University
% michael.rabbat@mcgill.ca
% 21 June 2007

ll = 0;
T = length(X);
for m=1:T
	Nm = length(X{m});
	p = 0;
	if (Nm <= 10)
		% Explicit marginalization over all permutations
		for t=1:Nm
			tau = zeros(Nm,1);
			p = permrecurse(t,1,Nm,tau,X{m},p,A,Pi);
		end
	else
		% Monte Carlo approximation
		% TODO: Fill this in later
		error('Monte Carlo for long paths has not been implemented yet');
	end
	if (p == 0)
		% Hmm, something's wrong here
		error(['Path ' num2str(m) ' has zero likelihood for this A and Pi']);
	end
	ll = ll + log(p) - log(gamma(Nm+1));
end

return;


% Function to recursively evaluate all permutations of a path
function p = permrecurse(visit, level, Nm, tau, x, p, A, Pi)

tau(visit) = level;
if (level == Nm)
	% Compute this likelihood term
	prob = Pi(x(tau(1)));
	for t=2:Nm
		prob = prob * A(x(tau(t-1)),x(tau(t)));
	end
	
	p = p + prob;
else
	% Continue recursive DF traversal of permutations
	for t=1:Nm
		if (tau(t)==0)
			p = permrecurse(t,level+1,Nm,tau,x,p,A,Pi);
		end
	end
end
tau(visit) = 0;

return;
