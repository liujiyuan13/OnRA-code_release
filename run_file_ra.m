clear
clc
warning off

proj_path = 'D:\Work\OnRA-code_release\';
addpath([proj_path, '/eval']);

data_path = 'D:\Work\datasets\mData\Fmatrix\';
data_names = {'Flower17', 'Caltech101-all', 'AwA', 'MNIST', 'YtVideo_sel'};

for id = 1
    
    data_name = data_names{id};
    fprintf('\n# data_name: %s', data_name);
    load([data_path, data_name, '_fea'], 'X', 'Y');
    k = length(unique(Y));
    V = size(X, 1);
    n = size(X{1}, 2);

    % compute kernels
    % normalize data
    for v=1:V
        X{v} = zscore(X{v}');
    end

    % parameters
    lambda_set = 10.^[-5:1:5];
    lambda_set = [0, lambda_set];
    n_iters = 10;
    maxiter = 100;

    metrics_meaning = {'acc'; 'nmi'; 'purity'; 'AR'; 'RI'; 'MI'; 'HI'; 'fscore'; 'precision'; 'recall'};
    
    res = zeros(length(lambda_set), n_iters, length(metrics_meaning));
    runtime = zeros(length(lambda_set), n_iters);
    objs = cell(length(lambda_set), n_iters);
    for i=1:length(lambda_set)
        parfor iter=1:n_iters
            % main algorithm
            tic;
            [H, obj] = mc_ra(X, k, lambda_set(i), maxiter);
            % final H 
            H = reshape(H, size(H,1), numel(H)/size(H,1));
            [y] = my_kmeans(H, k);
            ts = toc;
            [eval] = my_eval_y(y, Y);
            fprintf('\nlambda: %d, iter: %d, time: %f, acc: %f, nmi: %f, pur: %f', i, iter, ts, eval(1), eval(2), eval(3)); 
            res(i,iter,:) = eval;
            runtime(i, iter) = ts;
            objs{i, iter} = obj;
        end
    end
    
    save([proj_path, 'res_ra/', data_name, '_res.mat'], 'data_name', 'lambda_set', 'res', 'runtime', 'objs');

end
