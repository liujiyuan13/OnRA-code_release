function [H, obj_cell] = mc_ra(X, k, lbd, maxiter)
%MC_RA Multi-view Clustering with Representation Alignment

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


flag = 1; t = 1;
while(flag)
    % alternate block
    for v = 1:V
        % cal M
        XvTHv = X{v}'*H(:,:,v); XvXvTHv = X{v}*XvTHv;
        Hplusv = sum(H, 3) - H(:,:,v);
        M = 2*c1(v)*XvXvTHv + 2*(lbd/(V*k))*Hplusv;
        % update H
        [U, ~, VV] = svds(M, k, 'la');
        H(:,:,v) = U*VV';
    end
    % cal obj
    obj1(t) = 0; obj2(t) = 0;
    for v = 1:V
        XvTHv = X{v}'*H(:,:,v);
        obj1(t) = obj1(t) + (1/V) * trace(XvTHv'*XvTHv) * c1(v);
        Hsum = sum(H, 3) - H(:,:,v);
        obj2(t) = obj2(t) + (1/(V^2*k)) * trace(Hsum'*H(:,:,v));
    end
    obj = obj1 + lbd*obj2; obj_cell = {obj1, obj2};
    % check if convergent
    if t>1 && (abs((obj(t-1)-obj(t))/obj(t))<1e-5 || t>maxiter)
        flag = 0;
    end
    t = t+1;
end 


end

