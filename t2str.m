function str = t2str(time)
  type = 's';
  
  if(time/60>1)
    time =  time/60;
    type = 'm';
    
    if(time/60>1)
      time =  time/60;
      type = 'h';  
    end 
  end

str = [num2str(round(time*1000)/1000) type];
end