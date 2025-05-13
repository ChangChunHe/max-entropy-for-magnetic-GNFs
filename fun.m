function res = fun(x)

global coef n_ele 

new_coef = coef;

x = x/n_ele;
idx = find(x>0);
x = x(idx);
new_coef = new_coef(idx);

res = sum(new_coef.*x.*log(x));

% res = sum(new_coef.*log(x.^2)); 
