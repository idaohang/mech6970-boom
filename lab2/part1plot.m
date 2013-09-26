function out=part1plot(data,svnum)
plotdatax=[];
plotdatay=[];
plotdataz=[];
for j=1:size(data,1)
   if data{j,1}==svnum
       plotdatax(end+1)=data{j,4}(2);
       plotdatay(end+1)=data{j,4}(1);
       plotdataz(end+1)=data{j,4}(3);
   end
end


scatter(plotdatax,plotdatay,20,[rand,rand,rand],'fill')