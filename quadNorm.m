function n = quadNorm(coordinates, elements,varargin)
%
% norm = quadNorm(coordinates, elements)
% norm = quadNorm(coordinates, elements, 'FLAG')
%
% Diese Funktion Berechnet die Orthogonalen mit Laenge 1 über alle Flächen
% FLAG:
% w -> Laenge entspricht Flaecheninhalt
%
% P.Schaefer

%% Parameterueberpruefung
w = 1;

optargin = size(varargin,2);
if(optargin>1)
    error('Zu viele Argumente');
elseif(optargin==1)
    if(ismember('w',varargin{1}))
        w = 0;
    end
end

    %Numbers needed
    s_ele = size(elements,1);
    

%%  calculate the Norm
    n = zeros(size(elements,1),3);
    for i = 1:s_ele
        % normalized Vector on every triangle
        tri = elements(i,:);
        a = (coordinates(tri(2),:)-coordinates(tri(1),:));
        b = (coordinates(tri(4),:)-coordinates(tri(1),:));
        N = cross(a',b');
        if(w)
            N = N/norm(N);
        end
        n(i,:) = N;
    end

end



