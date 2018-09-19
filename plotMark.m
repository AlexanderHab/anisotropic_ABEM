function  plotMark(ele,coordinates,elements,varargin)
% plotMark(ele,coordinates,elements)


type = 'xg';
if(size(varargin,2)==1)
    type = varargin{1};
end

hold on;
%   disp 'Plot Updated'
for idx = ele
  if(length(idx)==1)
    current = sum(coordinates(elements(idx,[2,4])',:),1)/2;
    scatter3(current(1),current(2),current(3),type);
  else
    current= [sum(coordinates(elements(idx(1),[2,4])',:),1)/2; ...
      sum(coordinates(elements(idx(2),[2,4])',:),1)/2];
    plot3(current(:,1),current(:,2),current(:,3),'r');
  end
end
hold off;

end

