clear
close all

%% Example 1

load exmpl_2DQuad

%export_mesh(coordinates,elements,neigh,[],'exmpl11')


marked = [2 3 1 3];
[coordinates elements neigh f2s sites] =...
  refineQuad(coordinates,elements,neigh,sites,marked);

export_mesh(coordinates,elements,neigh,f2s,'exmpl12')

marked = ones(1,9);
marked(7) = 2;
% marked([4,7]) = 3;
[coordinates elements neigh f2s sites] =...
  refineQuad(coordinates,elements,neigh,sites,marked);

export_mesh(coordinates,elements,neigh,f2s,'exmpl13')

plotShape(coordinates,elements([4 10 13 12 9] ,:),'tb',...
  {'(4)' '(10)' '(13)' '(12)' '(9)'})
ylim([-.05 0.67])
xlim([-.05 0.55])
view(2)
text(.07,.27,0,'1');
text(.17,.27,0,'5');

text(.23,.33,0,'2');
text(.23,.43,0,'6');

text(.17,.47,0,'3');
text(.07,.47,0,'7');

text(.01,.43,0,'4');
text(.01,.33,0,'8');

print('-r600','-depsc',['../doc/fig/exmpl13_nei_part.eps'])

fid = fopen(['../doc/fig/exmpl13_nei_part.tex'],'w');
% fprintf(fid,'\\begin{figure}\n');
fprintf(fid,'\\begin{tabular}{>{\\columncolor{gray}}rcccccccc}\n');
fprintf(fid,'\\rowcolor{gray}\n Index & n1 & n2 & n3 & n4 & n5 & n6 & n7 & n8\\\\');

[m n] =size(neigh);

i = 4;
fprintf(fid,'\n   %d & %d & %d & %d & %d & %d & %d & %d & %d',i,neigh(i,1),neigh(i,2),neigh(i,3),neigh(i,4),neigh(i,5),neigh(i,6),neigh(i,7),neigh(i,8));

fprintf(fid,'\n\\end{tabular}');
% fprintf(fid,['\n\\caption{' file ' - Nachbarn}']);
fprintf(fid,['\n\\label{exmpl13:nei:part:nei}']);
% fprintf(fid,'\n\\end{figure}');
fclose(fid);

%% Example Right&Wrong

plotShape(coordinates,elements(:,:),'be');
view(2);
axis off;
print('-r600','-depsc',['../doc/fig/net1.eps'])
% system(['epstopdf ../doc/fig/net1.eps']);

load exmpl_2DQuad_wrong

plotShape(coordinates,elements(:,:),'be');
view(2);
axis off;
print('-r600','-depsc',['../doc/fig/net_wrong.eps'])
% system(['epstopdf ../doc/fig/net_wrong.eps']);

%% Example Single Element
load exmpl_2DQuad2

plotShape(coordinates,elements(:,:),'ben');
hold on
text(-.05,-.05,0,'k1');
text(1.05,0,0,'k2');
text(1.05,1.05,0,'k3');
text(-.05,1.05,0,'k4');

text(.5,-.05,0,'e1');
text(1.05,.5,0,'e2');
text(.55,1.1,0,'e3');
text(-.05,.5,0,'e4');

text(.55,.5,.7,'n');
hold off


% view(2);
% axis off;
print('-r600','-depsc',['../doc/fig/net_single.eps'])
% system(['epstopdf ../doc/fig/net_single.eps']);

%% Example Objects
load exmpl_2DLShape
export_mesh(coordinates,elements,neigh,[],'exmpl_2DLShape')

load exmpl_2DQuad
export_mesh(coordinates,elements,neigh,[],'exmpl_2DQuad')

load meshSave/1t05n05_2DQuad_12
plotShape(coordinates,elements,'b');
view(2)
print('-r600','-depsc','../doc/fig/exmpl_2DQuad_12')

load exmpl_3DCube
export_mesh(coordinates,elements,neigh,[],'exmpl_3DCube')

load exmpl_3DFichCube
export_mesh(coordinates,elements,neigh,[],'exmpl_3DFichCube')

load meshSave/1t05n05_3DFichCube_6
plotShape(coordinates,elements,'bv');
print('-r600','-depsc','../doc/fig/exmpl_3DFichCube_6')


%% Example refType
file = 'refType';
input = 'exmpl_2DQuad2';
load(input);

plotShape(coordinates,elements(:,:),'bet',{'T'});

view(2);
axis off;
% set(gcf, 'paperpositionmode', 'auto');
print('-r600','-depsc',['../doc/fig/' file '_full.eps'])
% system(['epstopdf ../doc/fig/' file '_full.eps']);


for i= 1:4
  
load(input)
figure(i)
[coordinates elements neigh f2s sites] = ...
  refineQuad(coordinates,elements,neigh,sites,i);

plotShape(coordinates,elements,'bet',{'T1','T2','T3','T4'});
view(2);
% axis off;
% set(gcf, 'paperpositionmode', 'auto');
print('-r600','-depsc',['../doc/fig/' file '_' num2str(i) '.eps'])
% system(['epstopdf ../doc/fig/' file '_' num2str(i) '.eps']);

end

%% voll Analytisch
% A_plots({'meshSave/1t05n05_3DFichCube_23'},'../doc/fig/1t05n05_3DFichCube')
% A_plots({'meshSave/1t05n05_3DCube_25'},'../doc/fig/1t05n05_3DCube')
% A_plots({'meshSave/1t05n05_2DQuad_32'},'../doc/fig/1t05n05_2DQuad')
% A_plots({'meshSave/1t05n05_2DLShape_30'},'../doc/fig/1t05n05_2DLShape')

%% Isotrop Uniform
A_plots({'meshSave/1t1n0_2DQuad_6'},'../doc/fig/1t1n0_2DQuad')
A_plots({'meshSave/1t05n0_2DQuad_17'},'../doc/fig/1t05n0_2DQuad')

A_plots({'meshSave/1t1n0_2DQuad_6',...
  'meshSave/1t05n0_2DQuad_16',...
  'meshSave/1t05n05_2DQuad_29'},'../doc/fig/1tn_2DQuad')

%% Semianalytisch
A_plots({'meshSave/132t05n05_2DQuad_30'},'../doc/fig/132t05n05_2DQuad')
A_plots({'meshSave/142t05n05_2DQuad_30'},'../doc/fig/142t05n05_2DQuad')
A_plots({'meshSave/142t05n05_3DFichCube_22'},'../doc/fig/132t05n05_3DFichCube')

% A_plots({'meshSave/1432t05n05_3DFichCube_21'},'../doc/fig/1432t05n05_3DFichCube')
% A_plots({'meshSave/1432t05n05_2DQuad_29'},'../doc/fig/1432t05n05_2DQuad')

% A_plots({'meshSave/2222t05n05_2DQuad_30'},'../doc/fig/2222t05n05_2DQuad')


close all