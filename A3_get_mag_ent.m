function [x,ent, uniq_bond] = A3_get_mag_ent(filename)

% filename = 'D:\work\research\magnetic-C\test_AFM\fig5\a\t4.vasp';

global coef n_ele

[C,S,atoms] = read_poscar(filename);
idx1 = str2num(atoms{2,1});
idx2 = str2num(atoms{2,2});
S = S(1:idx1,:);
S = S*C;
S(:,3) = 10;
% perm_table = get_perms(S,0.5);
perm_table = 1:size(S,1);
bond = [];
d = zeros(size(S,1));
for ii = 1:size(S,1)
    for jj = ii+1:size(S,1)
        d(ii,jj) = norm(S(ii,:) - S(jj,:));
        if d(ii,jj)<2.1
            bond = [bond;ii jj];
        end
    end
end
d = d + d';

for ii = 1:size(bond,1)
    t = bond(ii,1:2);
    min_seq = sortrows(sort(perm_table(:,t),2));
    bond(ii,3:4) = min_seq(1,:);
end

uniq_bond = unique(bond(:,3:4),'rows');
num_vars = size(uniq_bond,1);
A = zeros(size(S,1)+1,num_vars);

for ii = 1:size(S,1)
    [a,b] =  find(bond(:,1:2) == ii);
    for k = 1:length(a)
        [~,idx] = ismember(bond(a(k),3:4),uniq_bond,'rows');
        A(ii,idx) = A(ii,idx) + 1;
    end
end

for ii = 1:size(bond,1)
    [~,idx] = ismember(bond(ii,3:4),uniq_bond,'rows');
    A(end,idx) = A(end,idx) + 1;
end
b = zeros(size(S,1)+1,1);
for ii = 1:size(S,1)
    t = length(find(d(ii,:)>0 & d(ii,:)<1.7));
    if t == 2
        b(ii) = 3;
    elseif t == 3
        b(ii) = 4;
    end
end
idx_2 = length(find(b==3));
b(end) = (4*size(S,1)-idx_2)/2;
[~,idx] = unique(A,'rows');
A = A(idx,:);
b = b(idx);

[n_ele,idx] = max(b);
coef = A(idx,:);
num_vars = size(A,2);
x0 = rand(1,num_vars);

options = optimoptions('fmincon','Algorithm','sqp-legacy',...
    'MaxFunctionEvaluations',300000,'OptimalityTolerance',1e-7);

[x,ent] = fmincon(@fun,x0,[],[],A,b,ones(num_vars,1),2*ones(num_vars,1),[],options);

% res = [x' uniq_bond ];

