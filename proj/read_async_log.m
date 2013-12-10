clear all
close all
clc
% profile on
location = ['..' filesep 'data' filesep 'final_proj_data' filesep];
filename = 'run_north.alog';
matFile = ['run_north_alog.mat'];
tic

pb_count = 0;

%% Open File, and figure out number of lines

use_pbar = false;

if use_pbar
    fid = fopen(sprintf('%s%s',location,filename),'rt');
    nLines = 0;
    while (fgets(fid) ~= -1),
      nLines = nLines+1;
    end
    fclose(fid);

    update_progress_bar = round(linspace(1,nLines,10));
    
end

%% Open file and read data from it

fid = fopen(sprintf('%s%s',location,filename));

disp('Reading in data from file.')

%Delete initial comments in file, initialze D, sadly we lose the first line
%of the file... anyone wanna fix this?
D=cell(1,4);

while (isempty(D{1}))
    curLine=fgets(fid);
    D = textscan(curLine,'%n%s%s%s','CommentStyle','%','bufsize',100000); 
end



%% Create variables
%search through all of the recieved devices and create new variables names
%for these devices.  Additionally, look at the message give by this device,
%and create a variable name for it.
global data;

disp('Scanning file for Devices, Messages, and Data.')
disp('Go get some coffee...')
%initilaizing variables
%length_time = length(time); %value used for progress bar
input = {0};           %this array will contain the names of the various devices and messages
dvc_count = 0;
%data will contain the actual values that correspond to the input array
data(1,1).val = NaN;
init.val = NaN;
% data = init.val*ones(100000,30); 
data(100000,30).val = NaN;
time_data(1,1).time = NaN;
time_data(3,10000).time = NaN;
%init_amount = length_time;

while ~feof(fid)
    
    if use_pbar
        pb_count = pb_count+1;
        %upate progress bar
         if( (pb_count == update_progress_bar) == 0 )
        
          else
            progressbar(pb_count/ nLines);
         end
    end
    

    
    curLine=fgets(fid);
    
    D = textscan(curLine,'%n%s%s%n','CommentStyle','%'); 
    time    = D{1};
    message = D{2}(:);
    device  = D{3}(:);
    val     = D{4}(:);
    
    %Remove this line
    value = val;
    
    if isempty(val) %val will be empty if an array was attempted to be read in
        %comment sytle removes [rxc] from the begining of matrixes
        %whitespace removes '{' and '}' from data
        D = textscan(curLine,'%n%s%s%s','CommentStyle',{'[',']'}, 'bufsize',100000,'Whitespace','}{ '); 
        time    =  D{1};
        message = D{2}(:);
        device  = D{3}(:);
        val     = D{4}(:);
    
        if isempty(val) %If val is still empty, then there must be nothing to read in...
            continue    %if so, just skip everything else, go on to the next file line...
        end
           
            

        %This converts val from a string to char to a matrix.  If this
        %suddnely starts causing errors, look in old code, it was done
        %differently there...(was read in as a string as was a time vampire
        try
            m=textscan(val{1}(1:end-1),'%n','Delimiter',',');
        catch
            %This typically occurs with stationID with the Novatel
            disp('the value read in is empty, and will not be logged')
            continue;
        end %ends try

        value = m{1}(:);

    end



    %Iterate through the devices we've acquired so far and see if one
    %already exists in the array and returns the index of any match.  If no
    %match occurs, it will return 0
    dvc_idx = 0;
    for j=1:dvc_count+1
        if isempty(input{j,1})
            %if this index is empty put the new vlaue here
            break;
        else
            if length(input{j,1}) == length(device{1})
                if input{j,1} == device{1}
                    dvc_idx = j;
                end
            end
        end
    end



    if dvc_idx > 0 %check to see if device already exists
        
    

        %This will iterate through the messages acquired so far. The index
        %of any match will be found.  It will return 0 if no match is found
        msg_idx = 2;
        data_msg_idx = 0;
        for j=2:length(input(dvc_idx,:))
            if isempty(input{dvc_idx,j})
                %if the value is empty put the new value here
                msg_idx = j;
                break;           %no sense in contining to scan
            else
                if length(input{dvc_idx,j}) == length(message{1})
                    if input{dvc_idx,j} == message{1};
                       msg_idx = -1; 
                       data_msg_idx = j;
                       break;     %no sense in contining to scan
                    end
                end
            end
        end

        %if no match is found add it to the array
        if msg_idx > -1 || isempty(input(dvc_idx,2))
            %adds new value to message structure
            if ~isempty(input{dvc_idx,msg_idx})
                %cell not empty--&gt; expand length
                 input{dvc_idx,length(input(dvc_idx,:))+1} = message{1};
                 
                 %Close Up gaps in non-cell data
                 if length( value )==1
                     close_up_gaps(dvc_idx, length(input(dvc_idx,:)), length(time_data(dvc_idx,2).time));
                 end
                
                 if time_data(dvc_idx,2).time(end) == time(1)                      
                     try
                         data (dvc_idx,length(input(dvc_idx,:))).val( length(time_data(dvc_idx,2).time) ) = value; %double up so it will fail on try...this is dumb i know, but I'm also very lazy
                     catch
                         data (dvc_idx,length(input(dvc_idx,:))).val{ length(time_data(dvc_idx,2).time) } = value;
                     end
                 else
                     time_data(dvc_idx,2).time(end+1) = time(1);                     
                     try                         
                         data (dvc_idx,length(input(dvc_idx,:))).val( length(time_data(dvc_idx,2).time) ) = value; 
                     catch
                         data (dvc_idx,length(input(dvc_idx,:))).val{ length(time_data(dvc_idx,2).time) } = value; 
                     end
                 end
                 

            else
                %Close Up gaps in non-cell data
                if length( value )==1
                    close_up_gaps(dvc_idx, msg_idx, length(time_data(dvc_idx,2).time));
                end
                
                %cell empty --&gt; add message
                input{dvc_idx,msg_idx} = message{1};
%                 data (dvc_idx,msg_idx).val{1} = value;
%                 data (dvc_idx,msg_idx).time = time(1);
                 if time_data(dvc_idx,2).time(end) == time(1) 
                     try                         
                         data (dvc_idx,msg_idx).val( length(time_data(dvc_idx,2).time) ) = value;                         
                     catch                         
                         data (dvc_idx,msg_idx).val{ length(time_data(dvc_idx,2).time) } = value;
                     end
                 else
                     time_data(dvc_idx,2).time(end+1) = time(1);                     
                     try                         
                         data (dvc_idx,msg_idx).val( length(time_data(dvc_idx,2).time) ) = value;  
                     catch
                         data (dvc_idx,msg_idx).val{ length(time_data(dvc_idx,2).time) } = value;
                     end
                 end

            end
        else
             %Close Up gaps in non-cell data
%              if length( value )==1
%                  close_up_gaps      (dvc_idx, data_msg_idx, length(time_data(dvc_idx,2).time));
%              end
             
            %don't add message, but add data
%             data(dvc_idx,data_msg_idx).val{end+1} = value;
%             data(dvc_idx,data_msg_idx).time(end+1) = time(1);
             if time_data(dvc_idx,2).time(end) == time(1) 
                 try
                     data (dvc_idx,data_msg_idx).val( length(time_data(dvc_idx,2).time) ) = value;  
                 catch
                     data (dvc_idx,data_msg_idx).val{ length(time_data(dvc_idx,2).time) } = value; 
                 end
             else
                 time_data(dvc_idx,2).time(end+1) = time(1);                 
                 try
                     data (dvc_idx,data_msg_idx).val( length(time_data(dvc_idx,2).time) ) = value;  
                 catch
                     data (dvc_idx,data_msg_idx).val{ length(time_data(dvc_idx,2).time) } = value;
                 end
             end

        end

    else
        %First time device is loaded
        dvc_count = dvc_count + 1;
        input{end+1,1} = device{1};
        input{end,2} = message{1};

        [r,c] = size(input);
        try
            data(r,2).val(1) = value;
            %data(r,2).val(2:init_amount) = NaN;
        catch
            %data(r,2).val{1} = value;
        end
        
        time_data(r,2).time = time(1); %this will be the only index where time is kept up with...

    end
    
   
      
end

if use_pbar
    progressbar(1);
    close
end
%% Create strucutres to put data into

disp('Creating structures and populating them')

dvc_names=cell(1,length(input(:,1))-1);
for i=2:length(input(:,1)) %changed start form 2 to 1
    
    %get only the message indexes that contain data
    msg_idx = 0;
    for j=1:length(input(i,:))
        if ~isempty(input{i,j})
            msg_idx = msg_idx + 1;
        end
    end
    
    %Generate variable names
    % look for brackets in the name - messes up the eval statement
    % things written to the DB with uMS have the pc name in [].
    idxStrt=findstr(input{i,1},'[');
    if (~isempty(idxStrt))
       input{i,1}=input{i,1}(1:idxStrt-1);
    end
    
    % look for dashes in device name
    idxStrt=strfind(input{i},'-');
    if (~isempty(idxStrt))
       input{i}(idxStrt)='_';
    end
    
    dvc_name      = genvarname(input(i,1));
    msg_name_temp = genvarname(input(i,2:msg_idx))';
    dvc_names{1,i}     = dvc_name;
    %Get the maximum length time vector
    max_time_idx = 2;
    max_time_length = length(time_data(i,2).time);


    
    %add in the time vector
    msg_name = {ones(length(data(i,msg_idx).val(:)),1) 'time' time_data(i,max_time_idx).time(:)};
    for k=1:msg_idx-1     %makes up for starting at i=2
        try
            msg_name = [msg_name msg_name_temp{k} {{data(i,k+1).val(1:max_time_length)'}}];
        catch
            msg_name = [msg_name msg_name_temp{k} {{data(i,k+1).val(:)}}];
        end

    end
    msg_name = [msg_name 'NaN' 'NaN' ];

    %actual creation of the data structures   
    curName = input{i,1};
    idxStrt=strfind(curName,'-');
    if (~isempty(idxStrt))
       curName(idxStrt)='_';
    end
    eval([curName '=struct(msg_name{2:end-1},NaN);'])
  
end

% clear name_strings
% leng_msg_name=length(msg_name);
% jj=0;
% name_strings=cell(1,floor(leng_msg_name/2)-1);
% for kk=2:2:leng_msg_name-2
%     jj=jj+1;
%     name_strings{1,(jj)}=msg_name{1,(kk)};
%     save(matFile,'name_strings','-append');
% end


for kk=2:length(dvc_names)
    if kk==2
        save(matFile, char(dvc_names{1,kk}));
    else
        save(matFile, char(dvc_names{1,kk}),'-append');
    end
end

% also save start and end times of the log file
save(matFile, 'time','-append');


%  profile viewer
%         profsave(profile('info'),'profile_results')
        
        
clear all used variables
 clear D c data data_msg_idx device dvc_count dvc_idx dvc_name ...
        fid filename i input j k message msg_idx msg_name ...
        msg_name_temp r time value ll max_time_idx max_time_length ...
        m val length_time dvc_names init kk matFile time_data ans ...
        init_amount


    
%print data structuers to the screen so that the user can see what was
%created
 whos
 toc

 
return;
