% Example parser
% Takes the address of a text file and parses the data
%
% The file should have a new line for each input with format:
% [input variable name][space][input value]
input_path = 'C:\Users\jimcc\Desktop\testData.txt';

% The following are the inputs in the code
% These assigned values serve as defaults
ip1 = 0;
ip2 = 0;
ip3 = 0;
ip4 = 0;

file = fopen(input_path,'r'); % Open file read only

while true
    thisLine = fgetl(file); % Read line
    if ~ischar(thisLine); break; end % Get out of here if at end
    words = strsplit(thisLine); % Split the input and value
    
    % Theres a problem if there weren't two values after split
    if size(words) ~= 2
        print('Input error!')
        break; 
    end
    
    % Now check to see which input the read one matches and assign
    % if an input is not a number the 'str2double' should be changed
    if strcmp(words(1), 'ip1')
        ip1 = str2double(words(2));
    end
    if strcmp(words(1), 'ip2')
        ip2 = str2double(words(2));
    end
    if strcmp(words(1), 'ip3')
        ip3 = str2double(words(2));
    end
end

print('Here are the new inputs')
ip1
ip2
ip3
ip4