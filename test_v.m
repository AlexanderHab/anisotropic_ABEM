
mex mex_build_V.cpp slpRectangle.cpp

%lade Netz
    load meshSave/test1_1_1

%Berechne V 2 mal
    %Voll Analytisch
%     V1 = mex_build_V(coordinates,elements([10,6],:),[0.7],[1]);
    %SemiAnalytisch
      V1 = mex_build_V(coordinates,elements([10,6],:),[0.5 ],[0]);
      V2 = build_A2(coordinates,elements([10,6],:));


%Vergleiche beide Matritzen
   Vdif = abs(V1-V2);
   figure(1)
   spy(Vdif>10^-16)
   max(max(Vdif))
   
%    figure(2)
%    max(max(abs(V2'-V2)))
%    spy(abs(V2'-V2)>10^-16)
   