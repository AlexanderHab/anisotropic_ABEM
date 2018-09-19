function [area, vector, sites] = areaQuad(coordinates, elements,varargin)
%
% [area, vector, sites] = area(coordinates, elements)
%
% Diese Funktion Berechnet den Flaecheninhalt zu jedem Element (area) und 
% zu jedem Element die Orthogonalen mit Laenge 1, des Weiteren werden auch
% die Seitenlaengen zu jedem Element gespeichert
%
% P.Schaefer

    vector = zeros(size(elements,1),3);
    area = zeros(size(elements,1),1);
    sites = zeros(size(elements,1),2);
    
%%  Flaecheninhalt und Normalen berechnen

    for i = 1:size(elements,1)
        % normalized Vector on every triangle
        tri = elements(i,:);
        a = (coordinates(tri(2),:)-coordinates(tri(1),:));
        b = (coordinates(tri(4),:)-coordinates(tri(1),:));
        N = cross(a',b');
        
%         area(i) = norm(N);
        
        N = N/norm(N);
        vector(i,:) = N;
        area(i) = norm(a) * norm(b);
        sites(i,:) = [norm(a) norm(b)];
        
    end

end



