function A_loadMesh(coordinates,elements,neigh)
% Laed ein Mesh
%
% A_loadMesh(coordinates,elements)

  global G_C;
  global G_E;
  global G_N;
  global G_T;
  global G_D;
  
  G_C = coordinates;
  G_E = elements;
  G_N = neigh;
  G_T = [];
  G_D = [];

end
