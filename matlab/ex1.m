% ex1.m - Small 2-by-2 network with measurements of length 4
%
% Michael Rabbat
% McGill University
% michael.rabbat@mcgill.ca
%
% 19 June 2007

clear;

% Parameters
n = 2; % number of nodes in the network
T = 100; % number of paths
Nm = 4; % number of nodes per path

%rand('state', 20070619);

A = rand(n);
for i=1:n
    A(i,:) = A(i,:)./sum(A(i,:));
end

Pi = rand(n,1);
Pi = Pi./sum(Pi);

% Generate some paths according to this Markov model
X = cell(T,1);
for m=1:T
	% Allocate space for the path
	X{m} = zeros(Nm,1);
	
	% Sample the starting node from Pi
	cumprobs = cumsum(Pi);
	larger = find(cumprobs >= rand);
	X{m}(1) = larger(1);
	
	% Sample remaining nodes in the path by taking a random walk
	for t=2:Nm
		if (sum(A(X{m}(t-1),:)) > 0)
			cumprobs = cumsum(A(X{m}(t-1),:));
			larger = find(cumprobs >= rand);
			X{m}(t) = larger(1);
		else
			% Reached a dead end, terminate the path here
			X{m} = X{m}(1:t-1);
			break;
		end
	end
end

% Shuffle these observations
Y = cell(T,1);
for m=1:T
	% Figure out the length of this path
	Nm = length(X{m});
	% Generate a random permutation
	[sorted tau] = sort(rand(Nm,1));
	% Shuffle
	Y{m} = X{m}(tau);
end

numTrials = 50;
API = cell(numTrials,2);
LogLik = zeros(numTrials,1);
L1err = zeros(numTrials,1);
for k=1:numTrials
	disp(['k = ' num2str(k)]);

	% Estimate the Markov chain parameters from shuffled data with no priors
	APri = zeros(n);
	PiPri = zeros(n,1);
	[Ahat,pihat] = nico(Y,n,APri,PiPri);

	API{k,1} = pihat;
	API{k,2} = Ahat;
	LogLik(k) = loglik(X,Ahat,pihat);
	L1err(k) = sum(sum(abs(Ahat - A))) + sum(abs(pihat - Pi));
end

figure(1), clf;
scatter(L1err,LogLik);
xlabel('$\ell_1$ Error', 'interpreter', 'latex');
ylabel('Log Likelihood');

