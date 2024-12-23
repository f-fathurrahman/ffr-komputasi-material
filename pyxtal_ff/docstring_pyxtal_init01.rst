PyXtal_FF develops Machine Learning Interatomic Potential.

Parameters
----------
descriptors: dict
    The atom-centered descriptors parameters are defined here.
    
    The list of the descriptors keys:
    - type: str (SO4)
        The type of atom-centered descriptors.
        + ACSF (BehlerParrinello Gaussian symmetry)
        + wACSF (weighted Gaussian symmetry)
        + EAD (embeded atom density)
        + SO4 (bispectrum)
        + SO3 (smoothed powerspectrum)
        + SNAP (similar to SO4 but the weighting and Rc schemes are adopted from LAMMPS)
    - Rc: float/dictionary
        The radial cutoff of the descriptors. Dictionary form is particularly for SNAP.
    - weights: dictionary
        The relative species weights.
    - N_train: int
        The number of crystal structures in training data set 
        to be converted into descriptors.
    - N_test: int
        The number of crystal structures in test data set 
        to be converted into descriptors.
    - ncpu: int
        The number of cpu core to use for converting crystal structures 
        into descriptors.
    - stress: bool (False)
        Compute rdxdr (needed for stress calculation) or not
    - force: bool (True)
        Compute dxdr (needed for force calculation) or not
    - cutoff: str
        The cutoff function.
    - parameters: dict
        Example,
        + BehlerParrinello
            {'G2': {'eta': [.3, 2.], 'Rs': [1., 2.]}, 
                'G4': {'eta': [.3, .7], 'lambda': [-1, 1], 'zeta': [.8, 2]}}
        + EAD
            {'L': 2, 'eta': [.3, .7], 'Rs': [.1, .2]}
        + SO4/Bispectrum
            {'lmax': 3}
        + SO3
            {'nmax': 1, 'lmax': 3, 'alpha': 2.0}
        + SNAP
            {'weights': {'Si': 1.0, 'O': 2.0},
                'Rc': {'Si': 4.0, 'O': 5.0}
    - zbl: dict
        {'inner': 4.0, 'outer': 4.5}

model: dict
    Machine learning parameters are defined here.

    The list of the model keys:
    - algorithm: str (*NN and *PR)
        The desired machine learning algorithm for potential 
        development. Choose between ['PolynomialRegression', 'PR'] or
        ['NeuralNetwork', 'NN'].
    - system: list of str (*NN, *PR, and *GPR)
        A list of atomic species in the crystal system.
        e.g. ['Si', 'O']
    - hiddenlayers: list of int (*NN)
        [3, 3] contains 2 hidden layers with 3 nodes each.
    - activation: str or list of str (*NN)
        The activation function for the neural network model.
        Currently, there are tanh, sigmoid, and linear.
    - random_seed: int (*NN)
        If the Neural Network is initialized from the beginning,
        the random_seed is used to established the initial 
        Neural Network weights randomly.
    - batch_size: int (*NN)
        batch_size is used for online learning. The weights of the 
        Neural Network is updated for each batch_size.
    - epoch: int (*NN and *GPR)
        A measure of the number of times all of the training vectors 
        are used once to update the weights.
    - device: str (*NN and *GPR)
        The device used to train: 'cpu' or 'cuda'.
    - force_coefficient: float (*NN, *PR, and *GPR)
        This parameter is used as the penalty parameter to scale 
        the force contribution relative to the energy.
    - stress_coefficient: float (*NN, *PR, and *GPR)
        This parameter is used as the balance parameter scaling
        the stress contribution relative to the energy.
    - stress_group: list of strings, not bool! (*NN, *PR, and *GPR)
        Only the intended group will be considered in stress training.
    - softmax_beta: float (*NN)
        The parameters for Softmax Energy Penalty function.
    - unit: str (*NN)
        The unit of energy ('eV' or 'Ha'). 
        The default unit of energy is 'eV'. If 'Ha' is used,
        Bohr is the unit length; otherwise, Angstrom is used.
    - restart: str (*NN)
        To continue Neural Network training from where it was left off.
    - optimizer: dict (*NN and *GPR)
        Define the optimization method used to update NN parameters.
    - path: str (*NN, *PR, *GPR)
        The user defined path to store the NN results.
        path has to be ended with '/'.
    - memory: str (*NN)
        There are two options: 'in' or 'out'. 'in' will use load all
        descriptors to memory as 'out' will call from disk as needed.
    - order: int (*PR)
        Order is used to determined the polynomial order. 
        For order = 1, linear is employed, and quadratic is employed 
        for order = 2.
    - d_max: int (*PR)
        The maximum number of descriptors used.
    - alpha: float (*NN and *PR)
        L2 penalty (regularization term) parameter.
    - norm: int (*PR)
        This argument defines a model to calculate the regularization
        term. It takes only 1 or 2 as its value: Manhattan or Euclidean
        norm, respectively. If alpha is None, norm is ignored.
    - noise: float (*GPR)
        The noise added to the Gaussian Kernel
    - kernel: str (*GPR)
        The kernel specifying the covariance function of the GPR.
        The current development allows "RBF" and "DotProduct".
        
(*) required.
(*NN) for Neural Network algorithm only.
(*PR) for Polynomial Regression algorithm only.
(*GPR) for Gaussian Process Regressor algorithm only.



 Command PyXtal_FF to run in 2 modes:
1. train
    In train mode, PyXtal_FF needs TrainData and/or TestData to be defined.
    
    TrainData: str
        TrainData indicates the location the training data set file.
        After the training of Neural Network is over, TrainData will be 
        self-evaluated for the accuracy.
    
    TestData: str
        TestData is an optional argument. If the test data set is present, 
        PyXtal_FF will evaluate the accuracy of the developed potential 
        with the TestData.


2. predict
    In predict mode, PyXtal_FF need the saved machine learning interatomic
    potential.

    mliap: str
        The machine learning interatomic potential.
