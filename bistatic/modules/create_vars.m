function create_vars(param)
    % Create variables in the caller workspace from struct fields of param
    
    if ~isstruct(param)
        error('Input must be a struct');
    end
    
    fieldsList = fieldnames(param);
    for k = 1:numel(fieldsList)
        varName = fieldsList{k};
        varValue = param.(varName);
        assignin('caller', varName, varValue);
    end
end