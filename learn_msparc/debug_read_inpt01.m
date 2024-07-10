clear all; close all; %#ok<*CLALL>
filename = 'TEMP_Si2_kpt_PBE0'; % prefix of input files

fid1 = fopen(strcat(filename, '.inpt'), 'r');
if( fid1 == -1 ) 
    error('\nCannot open file "%s.inpt"\n', filename);
end

ACE_FLAG = 0;

while(~feof(fid1))
    C_inpt = textscan(fid1, '%s', 1, 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
    str = string(C_inpt{:}); % for R2016b and later, can use string()
    % skip commented lines starting by '#'
    %if (isempty(str) || str(1) == '#' || strcmp(str,'undefined'))
    if (isempty(str) || str(1) == '#')
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
        fprintf('skipping current line!\n');
        continue;
    end

    if( strcmp(str, 'ACE_FLAG:') )
        disp(' ')
        disp('===== Begin reading ACE_FLAG')
        C_param = textscan(fid1, '%f', 1, 'delimiter', ' ', 'MultipleDelimsAsOne', 1);
        disp(['C_param = ', C_param]);
        disp(['Class C_param ', class(C_param)])
        ACE_FLAG = C_param{1};
        %S.ACEFlag = C_param{:}; % ffr
        disp(['*** Read ACE_FLAG = ', ACE_FLAG]);
        disp('===== End of reading ACE_FLAG')
        disp(' ')
        textscan(fid1, '%s', 1, 'delimiter', '\n', 'MultipleDelimsAsOne', 0); % skip current line
    end % if

end

