for i=1:10
   try 
       x=x+1;
       disp(num2str(x));
   catch
       x=1;
   end
end