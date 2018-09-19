function plotShape(coordinates, elements, varargin)
%
% plotShape(coordinates,elements)
% plotShape(coordinates,elements,'FLAG','VEC')
%
% Diese Funktion Zeichnet alle Vierecke mit ausgefuellten Flaechen.
% FLAG:
% c -> Koordinaten als rote Kreise darstellen
% b -> nur Kanten der Elemente einzeichnen
% n -> Normen auf ein Element einzeichnen (mit Laenge 1)
% a -> Normen auf ein Element einzeichnen (Laenge durch Flaecheninhalt)
% s -> Flaechen werden gefaerbt, wobei VEC die Farben der Elemente angibt
% t -> Element Text wird angezeigt (oder Element Index)
% d -> Koordinaten Text wird angezeigt (oder Koordinaten Index)
% e -> dicke linien fuer den Export
% v -> gestrichelte linien fuer unsichtbare
%
% P.Schaefer

%% Parameterueberpruefung
c = 0;
e = 1;
n = 0;
t = [0 0];
f = 0;
v = 0;
line = 1;
optargin = size(varargin,2);
if(optargin>2)
    error('Zu viele Argumente');
elseif(optargin>=1)
    if(ismember('b',varargin{1}))
        e = 0;
    elseif(ismember('s',varargin{1}))
        if(optargin==2 && length(varargin{2})==size(elements,1))
            e = 2;
            color = reshape(varargin{2},length(varargin{2}),1);
%             maxc = max(color);
%             minc = min(color);
%             color = (color - minc)/maxc;
%             col = [color ones(length(color),2) ]
            
        else
            error('plotShape:COLOR','Optionale Parameter fuer Faerbung sind Falsch')
        end
%     elseif(ismember('t',varargin{1}))
%         if(optargin==2 && length(varargin{2})==size(elements,1))
%             t = 1;
% %             color = reshape(varargin{2},length(varargin{2}),1);
% %             maxc = max(color);
% %             minc = min(color);
% %             color = (color - minc)/maxc;
% %             col = [color ones(length(color),2) ]
% %             
% %         else
% %             error('plotShape:COLOR','Optionale Parameter fuer Faerbung sind Falsch')
%         end
    end
    
    if(ismember('e',varargin{1}))
        line = 2;
    end    
    
    if(ismember('c',varargin{1}))
        c = 1;
    end
    
    if(ismember('v',varargin{1}))
        v = 1;
    end
    if(ismember('a',varargin{1}))
        n = 2;
    elseif(ismember('n',varargin{1}))
        n = 1;
    end
    if(ismember('t',varargin{1}))
        desc = {};
        t(1) = 1;
        if(optargin==2 && length(varargin{2})>=size(elements,1))
            desc = varargin{2};
        end
    end
    if(ismember('d',varargin{1}))
        desc = {};
        t(2) = 1;
    end
end

%% Koordinaten einzeichnen
if(c)
    scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'r');
    hold on
end

% [b eles] = unique(elements,'rows');
%% Flaechen
eles = 1:size(elements,1);
coos = 1:size(coordinates,1);

anorm = quadNorm(coordinates,elements);
anorm(:,1:2) = anorm(:,1:2)*-1;

if(e==1)
    for idx = eles
        current = coordinates(elements(idx,[1:4,1])',:);
%         current(3,:) = current(3,:)-current(1,:)+current(2,:);
        fill3(current(:,1),current(:,2),current(:,3),[0,204,102]/255); % Zeichnet Oberflaeche
        hold on
    end
elseif(e==2)
    for idx = eles
        current = coordinates(elements(idx,[1:4,1])',:);
%         [current(:,1),current(:,2),current(:,3)]
%         col(idx,:)
%         current(3,:) = current(3,:)-current(1,:)+current(2,:);
        fill3(current(:,1),current(:,2),current(:,3),color(idx)); % Zeichnet Oberflaeche
%         view(2)
        hold on
    end
else
    for idx = eles
        if(sum(anorm(idx,:)<0)&&v==1)
          ltype = ':';
        else
          ltype = '-';
        end
        current = coordinates(elements(idx,[1:4,1])',:);
%         current(3,:) = current(3,:)-current(1,:)+current(2,:);
        plot3(current(:,1),current(:,2),current(:,3),ltype,'LineWidth',line,'color',[0,102,102]/255); % Zeichnet nur Kanten
        hold on
    end
end

alpha(0.3)

%% Zusatz
if(n)
    if(n==2)
        anorm = quadNorm(coordinates,elements,'w');
    else
        anorm = quadNorm(coordinates,elements);
    end
    for idx = eles
        current = sum(coordinates(elements(idx,[2,4])',:),1)/2;
        %Zeichnet auch den ersten Punkt ein !
        current = [current ; current+anorm(idx,:)];
        if(f==1)
          current = [current ; coordinates(elements(idx,1)',:)];
        end
        
        plot3(current(:,1),current(:,2),current(:,3),'r'); % Zeichnet Oberflaeche
        scatter3(current(2,1),current(2,2),current(2,3),'xr');
        hold on
    end
end


anorm = quadNorm(coordinates,elements);
anorm(:,1:2) = anorm(:,1:2)*-1;


if(t(1))
    for idx = eles
        current = sum(coordinates(elements(idx,[2,4])',:),1)/2;
        if(sum(anorm(idx,:)<0)&&v==1)
            cola = [.6 .6 .6];
        else
            cola = [0 0 0];
        end
        if(isempty(desc))
            text(current(1),current(2),current(3),['(' num2str(idx) ')'],'color',cola);
        else
            text(current(1),current(2),current(3),[desc{idx}],'color',cola);
        end
        hold on
    end
 end

if(t(2))
    for idx = coos
        current = coordinates(idx,:);

        cola = 'r';
            
        if(isempty(desc))
            text(current(1)+0.01,current(2)+0.03,current(3)+0.01, num2str(idx) ,'color',cola);
        else
            text(current(1)+0.01,current(2)+0.03,current(3)+0.01,desc{idx},'color',cola);
        end
        hold on
    end
end


%     anorm = quadNorm(coordinates,elements);
%     for idx = 1:eles
%         current = sum(coordinates(elements(idx,[2,4])',:),1)/2;
%         current = [current ; current+anorm(idx,:)*ind(idx);coordinates(elements(idx,1)',:)];
%         plot3(current(:,1),current(:,2),current(:,3),'r'); % Zeichnet Oberflaeche
%         scatter3(current(2,1),current(2,2),current(2,3),'xr');
%         hold on
%     end

xdif = max(coordinates(:,1)) - min(coordinates(:,1));
ydif = max(coordinates(:,2)) - min(coordinates(:,2));
zdif = max(coordinates(:,3)) - min(coordinates(:,3));

if(ydif)
    ylim([min(coordinates(:,1))-ydif/10 max(coordinates(:,1))+ydif/10])
end
if(xdif)
    xlim([min(coordinates(:,2))-xdif/10 max(coordinates(:,2))+xdif/10])
end
if(zdif)
    zlim([min(coordinates(:,3))-zdif/10 max(coordinates(:,3))+zdif/10])
end
xlabel 'x'
ylabel 'y'
zlabel 'z'



hold off

end

