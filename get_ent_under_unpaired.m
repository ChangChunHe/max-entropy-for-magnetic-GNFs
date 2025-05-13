function y = get_ent_under_unpaired(filename,unpaired_ele)
% filename = 'D:\work\research\magnetic-C\test_AFM\fig5\a\t2.vasp';
coef_set = 0;

global coef n_ele %filename

[C,S,atoms] = read_poscar(filename);
idx1 = str2num(atoms{2,1});
S = S(1:idx1,:);
S = S*C;
% perm_table = get_perms(S,1);
perm_table = 1:size(S,1);
% load matlab
max_d = 1.9;
bond = [];
d = zeros(size(S,1));
for ii = 1:size(S,1)
    for jj = ii+1:size(S,1)
        d(ii,jj) = norm(S(ii,:) - S(jj,:));
        if d(ii,jj) < max_d
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

iso_atoms_idx = [];
b = zeros(size(S,1)+1,1);
for ii = 1:size(S,1)
    t = length(find(d(ii,:)>0 & d(ii,:)<max_d));
    if t == 2
        b(ii) = 3;
        iso_atoms_idx = [iso_atoms_idx ii];
    elseif t == 3
        b(ii) = 4;
    end
end

iso_atoms_lu = ones(1,size(S,1));
iso_atoms_lu(iso_atoms_idx) = 1;

idx_2 = length(find(b==3));
b(end) = (4*size(S,1)-idx_2)/2;

[ii,uniq_atom,kk] = unique(A,'rows');
f = histc(kk,1:numel(uniq_atom)); % Frequency
A = ii;
b = b(uniq_atom);
add_iso = 2*eye(size(A,1));
add_iso(end,:) = f;
A = [A add_iso];
A(:,end) = [];

iso_atoms_lu = iso_atoms_lu(uniq_atom(1:end-1));

coef = A(end,:);coef(end-size(add_iso,1)+2:end) = coef_set;
n_ele = b(end);

A = [A;ones(1,size(A,2))];
A(end,1:size(bond,1)) = 0;
b = [b;unpaired_ele];

x0 = [rand(1,size(uniq_bond,1))+1  rand(1,size(A,2)-size(uniq_bond,1))];
% fminconOptions = optimset('Display', 'iter-detailed', 'Algorithm', 'sqp', ...
%     'AlwaysHonorConstraints', 'bounds');
fminconOptions = optimset( 'Algorithm', 'sqp', ...
    'AlwaysHonorConstraints', 'bounds','MaxFunEvals',1000000,'MaxIter',100000);
[x,y] = fmincon(@fun,x0,[],[],A,b,...
    [1*ones(1,size(uniq_bond,1)) 0*iso_atoms_lu],...
    [2*ones(1,size(uniq_bond,1)) 1*iso_atoms_lu],...
    [],fminconOptions);

% err = sum(abs(A*x'-b));
% mag = sum(x(size(bond,1)+1:end));
