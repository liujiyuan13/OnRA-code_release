clear
clc

%%
figure(1)
set(gcf,'Position',[300 500 1300 350]);
color_list = {'#0072BD', '#D95319', '#EDB120', '#77AC30', '#A2142F'};

% first gaussian distribution
subplot1 = subplot('position', [0.05, 0.2, 0.15, 0.47]);
% set(gca, 'FontSize', 20);

n=1000; 
mu=[4.25 4.25]; 
Sigma=[1.5,0;0,1.5]; 
r1=mvnrnd(mu,Sigma,n);
scatter(r1(:,1),r1(:,2), '.');
hold on

% second gaussian distribution
n=1000; 
mu=[8 8]; 
Sigma=[1.5,0;0,1.5]; 
r2=mvnrnd(mu,Sigma,n);
scatter(r2(:,1),r2(:,2), '.');
hold on

xlim([0 12]);
ylim([0 12]);
% set(gca,'ytick',[0:500:2000],'yticklabel',[0:5:20]);
ylabel('2-nd Dimension', 'Fontsize', 10);
xlabel('1-st Dimension', 'Fontsize', 10);
title('Gaussian Distribution', 'Fontsize', 10.5);
% axis square
grid on

for iter = 1:10000
    % get ind to compute
    num_a = 100;
    ind1 = randperm(n);
    ind1 = ind1(1:num_a);
    
    ind2 = randperm(n);
    ind2 = ind2(1:num_a);
    
    r1_a1 = r1(ind1,:);
    r2_a1 = r2(ind1,:);
    
    r1_a2 = r1(ind2,:);
    r2_a2 = r2(ind2,:);
    
    % compute representation alignment
    rep_alignment(1,iter) = trace(r1_a1*r1_a2')/sqrt(trace(r1_a1*r1_a1')*trace(r1_a2*r1_a2'));
    rep_alignment(2,iter) = trace(r2_a1*r2_a2')/sqrt(trace(r2_a1*r2_a1')*trace(r2_a2*r2_a2'));
    rep_alignment(3,iter) = trace(r1_a1*r2_a1')/sqrt(trace(r1_a1*r1_a1')*trace(r2_a1*r2_a1'));
    
    % compute vector alignment
    tmp1 = diag(r1_a1*r1_a2')./sqrt(diag(r1_a1*r1_a1').*diag(r1_a2*r1_a2'));
    tmp2 = diag(r2_a1*r2_a2')./sqrt(diag(r2_a1*r2_a1').*diag(r2_a2*r2_a2'));
    tmp3 = diag(r1_a1*r2_a1')./sqrt(diag(r1_a1*r1_a1').*diag(r2_a1*r2_a1'));
    vec_alignment(1,iter) = mean(tmp1);
    vec_alignment(2,iter) = mean(tmp2);
    vec_alignment(3,iter) = mean(tmp3);
    
end

subplot2 = subplot('position', [0.25, 0.2, 0.15, 0.47]);
% set(gca, 'FontSize', 10);

n_bin = 100;
minmin = min(min(rep_alignment));
maxmax = max(max(rep_alignment));
edges = (minmin:(maxmax-minmin)/n_bin:maxmax);
h1_rep = histogram(rep_alignment(1,:), edges); hold on
h2_rep = histogram(rep_alignment(2,:), edges); hold on
h3_rep = histogram(rep_alignment(3,:), edges); hold on
ylim([0 1500]);
set(gca,'ytick',[0:500:2000],'yticklabel',[0:5:20]);
ylabel('Percentage (%)', 'Fontsize', 10);
xlabel('Quantity', 'Fontsize', 10);
title('Representation Alignment', 'Fontsize', 10.5);
legend('Distribution 1', 'Distribution 2', 'Cross Distribution'); 
grid on

subplot3 = subplot('position', [0.45, 0.2, 0.15, 0.47]);
% set(gca, 'FontSize', 10);

n_bin = 100;
minmin = min(min(vec_alignment));
maxmax = max(max(vec_alignment));
edges = (minmin:(maxmax-minmin)/n_bin:maxmax);
h1_vec = histogram(vec_alignment(1,:), edges); hold on
h2_vec = histogram(vec_alignment(2,:), edges); hold on
h3_vec = histogram(vec_alignment(3,:), edges); hold on
ylim([0 1500]);
set(gca,'ytick',[0:500:2000],'yticklabel',[0:5:20]);
ylabel('Percentage (%)', 'Fontsize', 10);
xlabel('Quantity', 'Fontsize', 10);
title(' Average of Vector Inner-product', 'Fontsize', 10.5);
legend('Distribution 1', 'Distribution 2', 'Cross Distribution'); 
grid on




%%

% xs = [0:0.1:5];
% 
% for is = 1:length(xs)
% 
%     fprintf('\niter: %d', is);
%     
% x = xs(is);
%     
% % figure(1)
% % first gaussian distribution
% n=1000; 
% mu=[5 5]; 
% Sigma=[1,0;0,1]; 
% r1=mvnrnd(mu,Sigma,n);
% % scatter(r1(:,1),r1(:,2));
% % hold on
% 
% % second gaussian distribution
% n=1000; 
% mu=[10-x 10-x]; 
% Sigma=[1,0;0,1]; 
% r2=mvnrnd(mu,Sigma,n);
% % scatter(r2(:,1),r2(:,2));
% % hold on
% 
% for iter = 1:10000
%     % get ind to compute
%     num_a = 100;
%     ind1 = randperm(n);
%     ind1 = ind1(1:num_a);
%     
%     ind2 = randperm(n);
%     ind2 = ind2(1:num_a);
%     
%     r1_a1 = r1(ind1,:);
%     r2_a1 = r2(ind1,:);
%     
%     r1_a2 = r1(ind2,:);
%     r2_a2 = r2(ind2,:);
%     
%     % compute representation alignment
%     rep_alignment(1,iter) = trace(r1_a1*r2_a1')/sqrt(trace(r1_a1*r1_a1')*trace(r2_a1*r2_a1'));
%     rep_alignment(2,iter) = trace(r1_a1*r1_a2')/sqrt(trace(r1_a1*r1_a1')*trace(r1_a2*r1_a2'));
%     rep_alignment(3,iter) = trace(r2_a1*r2_a2')/sqrt(trace(r2_a1*r2_a1')*trace(r2_a2*r2_a2'));
%     
%     % compute vector alignment
%     tmp1 = diag(r1_a1*r2_a1')./sqrt(diag(r1_a1*r1_a1').*diag(r2_a1*r2_a1'));
%     tmp2 = diag(r1_a1*r1_a2')./sqrt(diag(r1_a1*r1_a1').*diag(r1_a2*r1_a2'));
%     tmp3 = diag(r2_a1*r2_a2')./sqrt(diag(r2_a1*r2_a1').*diag(r2_a2*r2_a2'));
%     vec_alignment(1,iter) = mean(tmp1);
%     vec_alignment(2,iter) = mean(tmp2);
%     vec_alignment(3,iter) = mean(tmp3);
%     
% end
% 
% 
% % hist
% n_bin = 100;
% minmin = min(min(rep_alignment));
% maxmax = max(max(rep_alignment));
% edges = (minmin:(maxmax-minmin)/n_bin:maxmax);
% % edges = (0.9994:0.000005:0.9999) - 0.0000025;
% h1_rep = hist(rep_alignment(1,:), edges);
% h2_rep = hist(rep_alignment(2,:), edges);
% h3_rep = hist(rep_alignment(3,:), edges);
% 
% n_bin = 100;
% minmin = min(min(vec_alignment));
% maxmax = max(max(vec_alignment));
% edges = (minmin:(maxmax-minmin)/n_bin:maxmax);
% % edges = (0.99965:0.000003:0.99995) - 0.0000015;
% h1_vec = hist(vec_alignment(1,:), edges);
% h2_vec = hist(vec_alignment(2,:), edges);
% h3_vec = hist(vec_alignment(3,:), edges);
% 
% overlap12_rep(is) = sum(min([h1_rep; h2_rep]));
% overlap13_rep(is) = sum(min([h1_rep; h3_rep]));
% overlap23_rep(is) = sum(min([h2_rep; h3_rep]));
% 
% overlap12_vec(is) = sum(min([h1_vec; h2_vec]));
% overlap13_vec(is) = sum(min([h1_vec; h3_vec]));
% overlap23_vec(is) = sum(min([h2_vec; h3_vec]));
% 
% end
% 
% figure(5)
% plot(xs, overlap13_rep); hold on
% plot(xs, overlap13_vec); hold on
% legend('rep', 'vec');

