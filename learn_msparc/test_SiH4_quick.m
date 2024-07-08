setup_path;

S = initialization('TEMP_SiH4_quick');
S.parallel = 0;  % no parallelization
[~,~,S] = electronicGroundStateAtomicForce(S.Atoms,S);

