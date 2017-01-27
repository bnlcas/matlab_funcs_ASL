function [collapsed_regress_matrix] = Collapse_Regression_Matrix(collapser, block_on_trigger, regress_matrix)
%% Divide regression matrix into smaller blocks, divided when the collapser input is triggered
% this could also be stimOnset triggered by {'S'} or trans/lex...

% block trigger is one if any of the blockontrigger tags appear
block_trigger = false(length(collapser),1); 
for i = 1:length(block_on_trigger)
    block_trigger = (block_trigger | strcmpi(collapser, block_on_trigger(i)));
end

% Get a listing of the Boolean Colums of the Regression Matrix:
% bool_cols is a list of the indecies of the boolean dummy variables in the regression matrix to be reduced
bool_cols = [];
j = 1;
for i = 1:size(regress_matrix,2)
    test = ((regress_matrix(:,i)==0) | (regress_matrix(:,i) == 1));
    if sum(~test) == 0
        bool_cols(j) = i;
        j = j+1;
    end
end

collapsed_regress_matrix = [];
i = 0;
j = 0;
while i<length(block_trigger)    
    i = i+1;
    if block_trigger(i)
        j = j+1;
        collapsed_regress_matrix(j,:) = regress_matrix(i,:);
        i = i+1;
        while ~block_trigger(i) & i < length(block_trigger)   
            collapsed_regress_matrix(j,bool_cols) = double(collapsed_regress_matrix(j,bool_cols) | regress_matrix(i,bool_cols));
            i = i+1;
        end
        i = i-1;
    end
end
