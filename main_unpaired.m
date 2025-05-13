clear
clc

unpaired_ele = linspace(0,1,101);
filename = '3.vasp';

for ii = 1:length(unpaired_ele)
    y(ii) = get_ent_under_unpaired(filename, unpaired_ele(ii));
end

y = y - y(1);
hold on
[~,idx] = min(y);
h = plot(unpaired_ele(idx), -y(idx),'o','MarkerSize',12);
c = get(h,'color');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(h,'markerfacecolor',c)

plot(unpaired_ele,-y,'LineWidth',2,'color',c)


box on
set(gca,'FontSize',18)
xlabel('$N_{\textrm{unpaired}}$')
ylabel('$S$')
