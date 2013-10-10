function [] = close_up_gaps( data_r, data_c, val_idx )

global data

try
   isempty( data(data_r,data_c) .val(val_idx-1) ) ;
catch
    try
        x_len = length( data(data_r,data_c).val );
        data(data_r,data_c).val(x_len+1:val_idx) = NaN; 
    catch
        data(data_r,data_c).val(1:val_idx) = NaN;  
    end
end
    
end

