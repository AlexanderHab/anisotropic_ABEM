function A_loadMeshF(file)
% Laed ein Mesh
%
% A_loadMesh(coordinates,elements)


  load(file)

  global G_C;
  global G_E;
  global G_N;
  global G_T;
  global G_D;
  global G_S;
  
  G_C = coordinates;
  G_E = elements;
  G_N = neigh;
  G_S = sites;
  G_T = [];
  G_D = [];
  
  if(exist('time','var')~=0)
    G_T = time;
  end
  if(exist('data','var')~=0)
    G_D = data;
  end
end



