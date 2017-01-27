function [filled_annot] = annot_fill(annot)
data_table = struct2table(annot);

%Very Poorly Written Code - only works when loops count from 1 up on
%indecies. Top must be filled in, since the function is implicitly
%recursive



for j = 1:size(data_table,2)
    col = table2array(data_table(:,j));   
    if isnumeric(col(1))
        isnull = (col == 0);   %zeros will be replaced
        for i = 2:size(data_table,1) % 2 necessary initialization condition
            if isnull(i)
                col(i) = col(i-1);
            end
        end
        data_table(:,j) = array2table(col);
    else
        isblank = (strcmpi(col,'') | strcmpi(col,' '));
        for i = 2:size(data_table,1)
            if isblank(i)
                col(i) = col(i-1);
            end
        end
        data_table(:,j) = col;
    end    
end

filled_annot = table2struct(data_table, 'ToScalar', true);

end