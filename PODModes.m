function k = PODModes(singular_values)
    % Change as desired, personally I prefer 100 - 1e-6
    threshold = 100 - 1.0e-6;
    
    % Calculate the Relative Info Content
    info_sum = sum(singular_values);
    info_cumsum = cumsum(singular_values);
    Relative_Information_Content = (info_cumsum / info_sum) * 100;

    % Choose the first k that satisfies the threshold
    k = find(Relative_Information_Content >= threshold, 1);
end