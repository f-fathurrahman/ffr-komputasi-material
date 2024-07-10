%
% Assumption: setup_path must be called first
%
clear all; close all; %#ok<*CLALL>

filename = 'TEMP_SiH4_quick'; % prefix of input files
% S = initialization();

% Set up inpt defaults
S = my_inpt_defaults();

% Read .inpt file
S = read_inpt(S, filename);

% Read .ion file
S = read_ion(S, filename);

% XXX Are there any new fields after calls to read_inpt and read_ion ?

% Read pseudopotential files
for ityp = 1:S.n_typ
    if (S.Atm(ityp).psptyp == 0)
        % WARNING: for TM format, lloc has to be set before reading
        S = ReadPseudoPot(S.Atm(ityp).lloc,S.Atm(ityp).typ);
        % ffr: This is not supported anymore??
        %
    elseif (S.Atm(ityp).psptyp == 1)
        S = readPseudopot(S, ityp, S.Atm(ityp).psdfname, S.Atm(ityp).typ);
    else
        error('Cannot recognize pseudopotential format!');
    end
end

% Set up more defaults based on input files
S = my_setup_defaults(S, filename);

% Calculate rb
S = my_calculate_rb(S);

% write initialized parameters into output file
S = my_write_output_init(S, filename);
% somehow output file will be referenced later

% This is needed
S.parallel = 0;  % no parallelization

