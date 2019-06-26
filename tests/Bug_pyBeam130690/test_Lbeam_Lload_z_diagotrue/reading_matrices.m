% Lettura Matrice di rigidezza (e M?)
filename = [pwd,'/pybeam_static_forstiffness.pch'];
formatSpec = '%1s%4s%7s%12s%16s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
wingpomatrix = cell2mat(raw);
clearvars filename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
keyboard
% gradi di liberta' indipendenti (??)
dim = 126;%ac.dof_ind;
K1 = sparse(dim,dim);
M1 = K1;

[m,~] = size(wingpomatrix);
row=0;
i=2;
column=0;
keyboard
while row<dim+1 && i<678   % subito dopo che finisce
    if isnan(wingpomatrix(i,4))
       row = (wingpomatrix(i,5)-1)*6+wingpomatrix(i,6);
    else 
        column = (wingpomatrix(i,4)-1)*6+wingpomatrix(i,5);
        K1(row,column) = wingpomatrix(i,6);
    end
    i=i+1;
    if column==dim
        i=i+1;
        break
    end
end
K = K1 + K1.';
K(1:size(K,1)+1:end) = diag(K1);
row = 0;
column=0;
%i=i+1;
keyboard
while row<dim+1 && i<2174923   % subito dopo che finisce
    if isnan(wingpomatrix(i,4))
       row = (wingpomatrix(i,5)-1)*6+wingpomatrix(i,6);
    else 
        column = (wingpomatrix(i,4)-1)*6+wingpomatrix(i,5);
        M1(row,column) = wingpomatrix(i,6);
    end
    i=i+1;
    if column==dim
        i=i+1;
        break
    end
end
M = M1 + M1.';
M(1:size(M,1)+1:end) = diag(M1);

 %clearvars -except K M
% save('wingPO_matrix')