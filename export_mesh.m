function export_mesh(coordinates, elements, neigh, f2s, file)

%% Daten
% save(file,'coordinates', 'elements','neigh','f2s')

%% Übersicht
plotShape(coordinates,elements,'bdtv');

if(length(unique(coordinates(:,3)))==1)
  view(2)
end
print('-r600','-depsc',['../doc/fig/' file '_ref.eps'])
% system(['epstopdf ../doc/fig/' file '_ref.eps']);

%% Koordinaten
fid = fopen(['../doc/fig/' file '_coo.tex'],'w');
% fprintf(fid,'\\begin{figure}\n');
fprintf(fid,'\\begin{tabular}{>{\\columncolor{gray}}rccc}\n');
fprintf(fid,'\\rowcolor{gray}\n Index & x1 & x2 & x3\\\\');

[m n] =size(coordinates);

for i = 1:m
    if(i~=1)
        fprintf(fid,'\\\\');
    end
    
    fprintf(fid,['\n    ' int2str(i) ...
      ' & ' num2str(coordinates(i,1),4)...
      ' & ' num2str(coordinates(i,2),4)...
      ' & ' num2str(coordinates(i,3),4)]);
end
fprintf(fid,'\n\\end{tabular}');
% fprintf(fid,['\n\\caption{' file ' - Koordinaten}']);
fprintf(fid,['\n\\label{' file ':coo}']);
% fprintf(fid,'\n\\end{figure}');
fprintf(fid,'\n');
fclose(fid);

%% Elemente
fid = fopen(['../doc/fig/' file '_ele.tex'],'w');
% fprintf(fid,'\\begin{figure}\n');
fprintf(fid,'\\begin{tabular}{>{\\columncolor{gray}}rcccc}\n');
fprintf(fid,'\\rowcolor{gray}\n Index & c1 & c2 & c3 & c4\\\\');

[m n] =size(elements);

for i = 1:m
    if(i~=1)
        fprintf(fid,'\\\\');
    end
    fprintf(fid,'\n   %d & %d & %d & %d & %d',i,elements(i,1),elements(i,2),elements(i,3),elements(i,4));
end
fprintf(fid,'\n\\end{tabular}');
% fprintf(fid,['\n\\caption{' file ' - Elemente}']);
fprintf(fid,['\n\\label{' file ':ele}']);
% fprintf(fid,'\n\\end{figure}');
fclose(fid);

%% Nachbarn
fid = fopen(['../doc/fig/' file '_nei.tex'],'w');
% fprintf(fid,'\\begin{figure}\n');
fprintf(fid,'\\begin{tabular}{>{\\columncolor{gray}}rcccccccc}\n');
fprintf(fid,'\\rowcolor{gray}\n Index & n1 & n2 & n3 & n4 & n5 & n6 & n7 & n8\\\\');

[m n] =size(neigh);

for i = 1:m
    if(i~=1)
        fprintf(fid,'\\\\');
    end
    fprintf(fid,'\n   %d & %d & %d & %d & %d & %d & %d & %d & %d',i,neigh(i,1),neigh(i,2),neigh(i,3),neigh(i,4),neigh(i,5),neigh(i,6),neigh(i,7),neigh(i,8));
end
fprintf(fid,'\n\\end{tabular}');
% fprintf(fid,['\n\\caption{' file ' - Nachbarn}']);
fprintf(fid,['\n\\label{' file ':nei}']);
% fprintf(fid,'\n\\end{figure}');
fclose(fid);

%% VaterSohn
fid = fopen(['../doc/fig/' file '_f2s.tex'],'w');
% fprintf(fid,'\\begin{figure}\n');
fprintf(fid,'\\begin{tabular}{>{\\columncolor{gray}}rcccc}\n');
fprintf(fid,'\\rowcolor{gray}\n Index & e1 & e2 & e3 & e4\\\\');

[m n] =size(f2s);

for i = 1:m
    if(i~=1)
        fprintf(fid,'\\\\');
    end
    fprintf(fid,'\n   %d & %d & %d & %d & %d',i,f2s(i,1),f2s(i,2),f2s(i,3),f2s(i,4));
end
fprintf(fid,'\n\\end{tabular}');
% fprintf(fid,['\n\\caption{' file ' - VaterSohn}']);
fprintf(fid,['\n\\label{' file ':f2s}']);
% fprintf(fid,'\n\\end{figure}');
fclose(fid);
end

