
file = 'refType';
input = 'exmpl_2DQuad2';
load(input);

plotShape(coordinates,elements(:,:),'bet',{'T'});

view(2);
axis off;
% set(gcf, 'paperpositionmode', 'auto');
print('-r600','-depsc',['../doc/fig/' file '_full.eps'])
system(['epstopdf ../doc/fig/' file '_full.eps']);


for i= 1:4
  
load(input)

[coordinates elements neigh f2s sites] = ...
  refineQuad(coordinates,elements,neigh,sites,i);

plotShape(coordinates,elements(:,:),'bet',{'T1','T2','T3','T4'});
view(2);
axis off;
% set(gcf, 'paperpositionmode', 'auto');
print('-r600','-depsc',['../doc/fig/' file '_' num2str(i) '.eps'])
system(['epstopdf ../doc/fig/' file '_' num2str(i) '.eps']);

end