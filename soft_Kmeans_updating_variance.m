%% Soft K-means algorithm for mixture of Gaussians
% Variance parameters and mixing coefficients are updated each 
% iteration using maximum likelihood results for univariate Gaussian.

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
SIG1 = [0.1 0; 0 22.9];
SIG2 = [0.1 0; 0 22.9];

coeffs= [0 0]; % mixture coefficients

%% Algorithm
% repeat until no changes in means
while (sum(sum(abs(means - old_means))) >1e-2)
    old_means = means;
    resp = zeros(N, 2); % vector of responsibilities
    
%% ASSIGNMENT STEP
% To increase k, add lines below (for K= 3-K) with form:
% 'gaussian(:,K) = mvnpdf([x1 x2], means(K,:), SIGK);

% generates gaussians for cluster k
        gaussian(:,1) = mvnpdf([x1 x2], means(1,:), SIG1);
        gaussian(:,2) = mvnpdf([x1 x2], means(2,:), SIG2);

    gaussianSum = sum(gaussian, 2); % sum of gaussians;

% This generates the responsibilities of each k for each data point n
    for n = 1:N
        for k = 1:K
         resp(n, k) = gaussian(n, k)./gaussianSum(n);
        end
    end
    
% This calculates the responsibilities: r_k(x-mu_k)^2 for each n
    diags= zeros(N, K);
    for n = 1:N
        for k = 1:K
            xn = data(n, :)-means(k, :);
            squared = xn'*xn; %Need to transpose matrix to square it
            diags(n, k) = squared(k,k); % only take the diagonal values
            topLine(n,k) = resp(n,k)*diags(n,k);
        end
    end
    
    top = sum(topLine, 1); %sum of all responsibilities*(x-mu)^2
    
% To increase k, add lines below (for K = 3-K) of the form:
% sigK = top/sum(resp(:,K))
% SIGK(1, 1) = sigK(1,1)
% SIGK(2, 2) = sigK(1,2)

    sig1 = top/sum(resp(:,1)); %ML for variance k=1
    sig2 = top/sum(resp(:,2)); %ML for variance k=2
    % place the variances in diagonal matrices (2x2)
    SIG1(1, 1) = sig1(1,1); 
    SIG1(2, 2) = sig1(1,2);
    SIG2(1, 1) = sig2(1,1);
    SIG2(2, 2) = sig2(1,2);
    
% update mixture coefficients.
    coeffs = (1/N)*sum(resp, 1);
 
    % Update means
    for k = 1:K
        for d = 1:D
            means(k,d) = sum(resp(:,k).*data(:,d));
        end
        means(k,:) = means(k,:)/sum(resp(:,k));
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
