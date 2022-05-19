function [H, obj_cell] = mc_cra(X, k, lbd, maxiter)
%MC_CRA Multi-view Clustering with Contrastive Representation Alignment

V = length(X);
n = size(X{1}, 1);

% Random init H
H = zeros(n, k, V);
Xm = cell2mat(X');
if size(Xm,2)<k
    Xm = [Xm, rand(n,k-size(Xm,2))];
end
[U, ~, ~] = svds(Xm, k, 'la');
for v = 1:V
    H(:,:,v) = U;
    XTX = X{v}'*X{v};
    c1(v) = 1 / sqrt(k * trace(XTX*XTX));
end
c2 = 1 / ((V-1)*k);

w_pos = 1/sqrt(2*n);
w_neg = 1/sqrt(2*n*(n-1));

flag = 1; t = 1;
while(flag)
    % alternate block
    for v = 1:V
        % cal M
        XvTHv = X{v}'*H(:,:,v); XvXvTHv = X{v}*XvTHv;
        BHv = (w_pos+w_neg)*H(:,:,v) - w_neg*repmat(sum(H(:,:,v),1), n, 1);
        GvHv = c1(v)*XvXvTHv;
        Hplusv = sum(H, 3) - H(:,:,v);
        BTHplusv = (w_pos+w_neg)*Hplusv - w_neg*repmat(sum(Hplusv,1), n, 1);
        M = 2*GvHv + 2*lbd*c2*BTHplusv;
        % update H
        [U, ~, VV] = svds(M, k, 'la');
        H(:,:,v) = U*VV';
    end
    % cal obj
    obj1(t) = 0; obj2(t) = 0;
    for v = 1:V
        XvTHv = X{v}'*H(:,:,v);
        obj1(t) = obj1(t) + (1/V) * trace(XvTHv'*XvTHv) * c1(v);
        Hplusv = sum(H, 3) - H(:,:,v);
        BTHplusv = (w_pos+w_neg)*Hplusv - w_neg*repmat(sum(Hplusv,1), n, 1);
        obj2(t) = (1/(V*(V-1)*k)) * trace(BTHplusv'*H(:,:,v));
    end
    obj = obj1 + lbd*obj2; obj_cell = {obj1, obj2};
    % check if convergent
    if t>1 && (abs((obj(t-1)-obj(t))/obj(t))<1e-5 || t>maxiter)
        flag = 0;
    end
    t = t+1;
end 


end

