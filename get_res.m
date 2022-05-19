clear
clc
warning off

proj_path = 'D:\Work\OnRA-code_release\';
addpath([proj_path, '/eval']);

data_path = 'D:\Work\datasets\mData\Fmatrix\';
data_names = {'Flower17', 'Caltech101-all', 'AwA', 'MNIST', 'YtVideo_sel'};

res_out_ra = zeros(length(data_names), 3);
res_out_cra = zeros(length(data_names), 3);

for id = 1:length(data_names)
    
    data_name = data_names{id};
    load([proj_path, 'res_ra/', data_name, '_res.mat'], 'res', 'objs', 'runtime');
    
    res = res * 100;
    res = squeeze(mean(res, 2));
    res_dac = res(2:end,:);
    res_dac = max(res_dac, [], 1);
    
    res_out_ra(id,:,1) = res(1,1:3);
    res_out_ra(id,:,2) = res_dac(1:3);
    
    load([proj_path, 'res_cra/', data_name, '_res.mat'], 'res', 'objs', 'runtime');
    
    res = res * 100;
    res_mean = squeeze(mean(res, 2));
    res_dac = res_mean(2:end,:);
    [res_dac, ind] = max(res_dac, [], 1);
    
    res_out_cra(id,:,1) = res_mean(1,1:3);
    res_out_cra(id,:,2) = res_dac(1:3);
    for i = 1:3
        std_out_cra(id,i,1) = std(res(1,:,i));
        std_out_cra(id,i,2) = std(res(1+ind(i),:,i));
    end
    
    res_out(id,:,1) = res_out_ra(id,:,1);
    res_out(id,:,2) = res_out_ra(id,:,2);
    res_out(id,:,3) = res_out_cra(id,:,2);
    
end