%% Soft K-means algorithms using full Gaussian mixture model.
% This is also known as the EM (Expectation-Maximization) Algorithm 
% for Gaussian Mixtures. Each Gaussian has a full covariance matrix,
% and each matrix has four parameters to be learned from the data

data = dlmread('old_faithful.dat', '',26, 1);

x1 = data( :, 1); % length of eruption
x2 = data( :, 2); % time since last eruption

K = 2; % code will need to be adjusted for different k values
N = length(data);
D = 2;

% randomly generate mean starting points
means =[(max(x1)-min(x1))*rand(K,1) +...% random starting means
    min(x1) (max(x2)-min(x2))*rand(K,1) + min(x2)];

old_means = zeros(K, D);

% starting point for means
SIG1 = cov(data);
SIG2 = SIG1;

% To increase K: make this a row vector with sum = 1.
% Each value can be 1/K for simplicity's sake
% Assign initial coefficients (randomly in this case)
PI_init = rand
coeff= [PI_init (1-PI_init)] % mixture coefficients need to sum to one

%% Algorithm
% repeat until no changes in means
while (sum(sum(abs(means - old_means))) >1e-2)
    old_means = means;
    resp = zeros(N, 2); % vector of responsibilities
    
%% ASSIGNMENT STEP
% To increase k, add lines below (for K= 3-K) with form:
% 'gaussian(:,K) = coeff(1,K)*(mvnpdf([x1 x2], means(K,:), SIGK));

% generates gaussians for cluster k
        gaussian(:,1) = coeff(1,1)*(mvnpdf([x1 x2], means(1,:), SIG1));
        gaussian(:,2) = coeff(1,2)*(mvnpdf([x1 x2], means(2,:), SIG2));

    gaussianSum = sum(gaussian, 2); % sum of gaussians;

% This generates the responsibilities of each k for each data point n
    for n = 1:N
        for k = 1:K
         resp(n, k) = gaussian(n, k)./gaussianSum(n);
        end
    end
    
    % for increased K, initialise more 2x2 'respSumK' matrices
    respSums1= [0,0 ; 0,0];
    respSums2= [0,0 ; 0,0];
    
% This calculates the SIGMA^2 for each k
    for n = 1:N
        for k = 1:K
            xn = data(n, :)-means(k, :);
            squared = xn'*xn; %Need to transpose matrix to square it
            t = resp(n,k)*squared;
            if k == 1
                respSums1 = respSums1 + t; % for increased K, add another
            elseif k == 2                  % if K == .... statement
                respSums2 = respSums2 + t; % and initalise another
            end                            %'respSumsK' 2x2 matrix 
    end
    
    %top = sum(t, 1); %sum of all responsibilities*(x-mu)^2
    
% To increase k, add lines below (for K = 3-K) of the form:
% SIGK = respSumsK/sum(resp(:,K))

    SIG1 = respSums1/sum(resp(:,1)); % for variance k=1
    SIG2 = respSums2/sum(resp(:,2)); % for variance k=2
    
% update mixture coefficients.
    coeffs = sum(resp, 1)/N;
 
    % Update means
    for k = 1:K
        for d = 1:D
            means(k,d) = sum(resp(:,k).*data(:,d));
        end
        means(k,:) = means(k,:)/sum(resp(:,k));
    end   
    end
end

%% Plot them
figure(1); clf
scatter(x1, x2); axis square, box on
xlabel('eruption duration(min)'); ylabel('time to next (min)')

figure(1); hold on
plot(means(:, 1), means(:,2), 'ksq', 'markersize', 15)

xx = linspace(0, 5);
yy = linspace(40, 100);
[XX, YY] = meshgrid(xx, yy);
 
f = mvnpdf([XX(:) YY(:)], means(1,:), SIG1);
f2 = mvnpdf([XX(:) YY(:)], means(2,:), SIG2);
contour(XX, YY, reshape(f,100,100))
contour(XX, YY, reshape(f2,100,100))
colorbar
