function recompute_time(filein,zeta)
%     dire = 'meshSave/';
    rows = 13;
    
%     filein = [dire '132t05n05_2DQuad_'];
    
%     zeta = {[2 2 2] [2 2 2] [2 2 2]};
   
for j=1:40
	disp(j);
    if(exist([filein num2str(j) '.mat'], 'file') ~= 2)
       break
    end
      

    load([filein num2str(j)]);
    [m n] = size(data);
    step = round(n/rows);
    
    typ = data(1,[2+(0:step-1)*rows])

    [coo_fine,ele_fine,neigh_fine,f2s,sit_fine]...
    =refineQuad(coordinates,elements,neigh,sites,2);

  for i = 1:length(typ)

      tic
      V_fine = mex_build_V(coo_fine,ele_fine,zeta{i},typ(i));
      time(j,i) = toc;
      clear V_fine
  end
%   time
%   data(1:j,[3+11+(0:step-1)*rows])
  data(1:j,[3+11+(0:step-1)*rows]) = time;
%   data(1:j,[3+11+(0:step-1)*rows])
  save([filein num2str(j)],'data','-append');
end

end
