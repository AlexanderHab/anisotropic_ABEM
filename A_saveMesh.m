function A_saveMesh(name)

global G_E;
global G_C;
global G_N;
global G_T;
global G_D;
global G_S;

time =[];
data =[];
coordinates = G_C;
elements = G_E;
neigh = G_N;
sites = G_S;
if(~isempty(G_T))
  time = G_T;
  data = G_D;
end


save (['meshSave/' name], 'coordinates', 'elements','neigh','time','data', 'sites')

end

