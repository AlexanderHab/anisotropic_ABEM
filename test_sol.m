
format longG

% Matrix MEX Funktion neu Compilieren
delete mex_build_V.mex*
%mex mex_build_V.cpp slpRectangle.cpp CFLAGS="\$CFLAGS -O2  -fopenmp" CXXFLAGS="\$CXXFLAGS -O2 -fopenmp" LDFLAGS="\$LDFLAGS -O2 -fopenmp"

%Alexx
mex mex_build_V.cpp slpRectangle.cpp CXXFLAGS="\$CXXFLAGS -O2 -fopenmp" 


% mex mex_build_V.cpp slpRectangle.cpp 

% Test ausführen

%Anzahl der Schritte oder wenn groeßer als 40 der Elemente
% steps = 10^4;

%Art der Berechnungen
 type = [8];
 zeta = {[1,1]};


%% Quad adaptiv anisotrop
%Datei
 %file = 'exmpl_3DFichCube';
 file = 'exmpl_3DCube';
% file = 'exmpl_2DQuad';
% file = 'meshSave/132t05n05_2DQuad_30';
 %file = 'meshSave/1432t05n05_3DFichCube_20';


steps =10^3; %Axx

%Adaptiv
 theta = 0.5;
%Anisotrop
 nu = 0.5;
 

tic
 [a, b ,fileo]=compute(file, steps, zeta, type, theta, nu, 0);
% 
 time = toc
 typeN = int2str(type);
 steps = size(a,1);

 A_plots({['meshSave/' fileo int2str(steps)]},...
     ['plots/' fileo ...
     int2str(steps)]);
  

%% Tests für export_example

% compute('exmpl_2DQuad', 10^4, {[2 2 2]}, [1], 1, 0, 0);
% compute('exmpl_2DQuad', 10^4, {[2 2 2]}, [1], .5, 0, 0);
% compute('exmpl_2DQuad', 29, {[2 2 2]}, [1], .5, .5, 0);

% compute('exmpl_3DFichCube', 1000, {[2 2 2]}, [1], .5, .5, 0);
% compute('exmpl_3DCube', 1000, {[2 2 2]}, [1], .5, .5, 0);
% compute('exmpl_2DLShape', 1000, {[2 2 2]}, [1], .5, .5, 0);

%compute('exmpl_2DQuad', 10^4, {[2 2 2] [2 2 2] [2 2 2]}, [1 4 2], .5, .5, 0);

% compute('exmpl_3DFichCube', 10^4, {[2 2 2] [2 2 2] [2 2 2]}, [1 4 2], .5, .5, 0);

% compute('exmpl_2DQuad', 10^4, {[0 2 2] [1 2 2] [2 2 2] [3 2 2]}, [2 2 2 2], .5, .5, 0);
