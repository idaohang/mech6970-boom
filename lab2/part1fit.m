function out=part1fit(data,svnum,destime)
    timedata=[];
    ecefx=[];
    ecefy=[];
    ecefz=[];
    for j=1:size(data,1)
        if data{1,1}=svnum
            timedata(end+1)=data{j,2};
            ecefx(end+1)=data{j,3}(1);
            ecefy(end+1)=data{j,3}(2);
            ecefz(end+1)=data{j,3}(3);

        end
    end
end