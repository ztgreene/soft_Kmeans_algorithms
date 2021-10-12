%% Soft K-means algorithm for two Guassians with simple fixed covariance matrix


data = dlmread('old_faithful.dat', '',26, 1);

x1 = data( :, 1); % length of eruption
x2 = data( :, 2); % time since last eruption

figure(1); clf
scatter(x1, x2); axis square, box on
xlabel('eruption duration(min)'); ylabel('time to next (min)')

MU1 = mean(x1); % means from sample
MU2 = mean(x2); %% don't need these

% Data has range of 3.5 in one dimension, and 53 in other
% (squared ratio of ranges is (53/3.5)^2 = 229
SIGMA = [0.1 0; 0 22.9]; % diagonal covariance matrix

K = 2;
N = length(data);
D = 2;

% randomly generate mean starting points
means =[(max(x1)-min(x1))*rand(K,1) +...% random starting means
    min(x1) (max(x2)-min(x2))*rand(K,1) + min(x2)];

old_means = zeros(K, D);

%% Algorithm
% repeat until no changes in means
while (sum(sum(abs(means - old_means))) >1e-2)
    old_means = means;
    resp = zeros(N, 2);
    topLine = zeros(N, 2);
% ASSIGNMENT STEP - update - change responsibility to Gaussian/sum of gaussians
% generates probs for cluster k
    for k = 1:K
        topLine(:,k) = mvnpdf([x1 x2], means(k,:), SIGMA);
    end

    bottomLine = sum(topLine, 2); % sum of gaussians;

% This generates the responsibilities of each k for each data point n
    for n = 1:N
        for k = 1:k
         resp(n, k) = topLine(n, k)./bottomLine(n);
        end
    end
    
    % Update
    for k = 1:K
        for d = 1:D
            means(k,d) = sum(resp(:,k).*data(:,d));
        end
        means(k,:) = means(k,:)/sum(resp(:,k));
    end
end

%% Plot them
figure(1); hold on
plot(means(:, 1), means(:,2), 'ksq', 'markersize', 15)

xx = linspace(0, 5);
yy = linspace(40, 100);
[XX, YY] = meshgrid(xx, yy);
 
f = mvnpdf([XX(:) YY(:)], means(1,:), SIGMA);
f2 = mvnpdf([XX(:) YY(:)], means(2,:), SIGMA);
contour(XX, YY, reshape(f,100,100))
contour(XX, YY, reshape(f2,100,100))
colorbar
