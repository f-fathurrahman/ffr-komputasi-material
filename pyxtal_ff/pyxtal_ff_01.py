import os
import torch

from utilities_database import Database
from models_polynomialregression import PR
from models_neuralnetwork import NeuralNetwork

__version__ = '0.2.3'

class PyXtal_FF():
    def __init__(self, descriptors=None, model=None, logo=True):
        if logo:
            self.print_logo()
        
        # Checking the keys in descriptors
        descriptors_keywords = ['type', 'Rc', 'weights', 'N_train', 'N_test', 'cutoff',
                                'force', 'stress', 'ncpu', 'parameters', 'base_potential']
        if descriptors is not None:
            for key in descriptors.keys():
                if key not in descriptors_keywords:
                    msg = f"Don't recognize {key} in descriptors. "+\
                          f"Here are the keywords: {descriptors_keywords}."
                    raise NotImplementedError(msg)

        # Set up default descriptors parameters
        self._descriptors = {'system': model['system'],
                             'type': 'Bispectrum',
                             'Rc': 5.0,
                             'weights': None,
                             'N': None,
                             'N_train': None,
                             'N_test': None,
                             'ncpu': 1,
                             'force': True,
                             'stress': True,
                             'cutoff': 'cosine',
                             'base_potential': False,
                             }
        
        # Update the default based on user-defined descriptors
        if descriptors is not None:
            self._descriptors.update(descriptors)
            if 'type' in descriptors and descriptors['type'] in ['EAD', 'ead']:
                _parameters = {'L': 3, 'eta': [0.1], 'Rs': [1.]}
            elif 'type' in descriptors and descriptors['type'] in ['SO3', 'SOAP']:
                _parameters = {'nmax': 3, 'lmax': 3, 'alpha': 2.0}
            else:
                _parameters = {'lmax': 3, 'rfac': 0.99363, 'normalize_U': False}
            if 'parameters' in descriptors:
                _parameters.update(descriptors['parameters'])
                self._descriptors['parameters'] = _parameters

            # Check for the SNAP type
            if self._descriptors['type'] in ['SNAP', 'snap']:
                if not isinstance(self._descriptors['weights'], dict):
                    msg = "The weights for SNAP type must be defined as a dictionary."
                    raise ValueError(msg)
                #if not isinstance(self._descriptors['Rc'], dict):
                #    msg = "The Rc for SNAP type must be defined as a dictionary."
                #    raise ValueError(msg)

        # Create new directory to dump all the results.
        # E.g. for default 'Si-O-Bispectrum/'
        if 'path' in model:
            self.path = model['path']
        else:
            _system = model['system']
            self.path = "-".join(_system) + "-"
            self.path += self._descriptors['type'] + "/"

        if logo:
            if not os.path.exists(self.path):
                os.mkdir(self.path)
            self.print_descriptors(self._descriptors)
        
        # Checking the keys in model.
        keywords = ['algorithm', 'system', 'hiddenlayers', 'activation', 
                    'random_seed', 'force_coefficient', 'unit', 'softmax_beta', 
                    'restart', 'optimizer', 'path', 'order', 'd_max', 
                    'epoch', 'device', 'alpha', 'batch_size', 'noise', 'kernel',
                    'norm', 'stress_coefficient', 'stress_group', 'memory']
        for key in model.keys():
            if key not in keywords:
                msg = f"Don't recognize {key} in model. "+\
                      f"Here are the keywords: {keywords}."
                raise NotImplementedError(msg)
        
        # Create model
        pr_keywords = ['PolynomialRegression', 'PR']
        nn_keywords = ['NeuralNetwork', 'NN']
        if 'algorithm' not in model:
            model['algorithm'] = 'NN'

        if model['algorithm'] in pr_keywords:
            self.algorithm = 'PR'
        elif model['algorithm'] in nn_keywords:
            self.algorithm = 'NN'
        else:
            msg = f"{model['algorithm']} is not implemented."
            raise NotImplementedError(msg)
        
        self._model = model

    def todict(self):
        return {"descriptor": self._descriptors, "model": self._model}

    def run(self, mode='train', TrainData=None, TestData=None, mliap=None):

        if mode == 'train':
            assert TrainData is not None, "TrainData can't be None for train mode."

            # Instantiate model
            print("Initializing model ...")
            self._MODEL(self._model)
            print("... Done initializing model")

            # Calculate descriptors.
            self._descriptors.update({'N': self._descriptors['N_train']})
            if not os.path.exists(self.path+'Train_db.dat') and not os.path.exists(self.path+'Train_db.db'):
                trainDB = Database(name=self.path+'Train_db')
                trainDB.store(TrainData, self._descriptors, True, self.path+'ase.db')
            else:
                trainDB = Database(name=self.path+'Train_db')
                trainDB.store(TrainData, self._descriptors, False)
            trainDB.close()

            if TestData is not None:
                EvaluateTest = True
                self._descriptors.update({'N': self._descriptors['N_test']}) 
                if not os.path.exists(self.path+'Test_db.dat'):
                    testDB = Database(name=self.path+'Test_db')
                    testDB.store(TestData, self._descriptors, True, self.path+'ase.db')
                else:
                    testDB = Database(name=self.path+'Test_db')
                    testDB.store(TestData, self._descriptors, False)
                testDB.close()

            else:
                EvaluateTest = False
            
            print("=========================== Training =============================\n")

            self.model.train('Train_db', optimizer=self.optimizer)
            self.model.save_checkpoint(des_info=self._descriptors)
            
            print("==================================================================\n")
            
            print(f"==================== Evaluating Training Set ====================\n")

            train_stat = self.model.evaluate('Train_db', figname='Train.png')
            
            print("==================================================================\n")

            if EvaluateTest:
                print("================= Evaluating Testing Set =====================\n")

                test_stat =  self.model.evaluate('Test_db', figname='Test.png')
                
                print("==============================================================\n")
            else:
                test_stat = None

            return (train_stat, test_stat)
        
        elif mode == 'predict':
            self._model['algorithm'] = torch.load(mliap)['algorithm']
            self.algorithm = self._model['algorithm']
            self._MODEL(self._model)
            self._descriptors = self.model.load_checkpoint(filename=mliap)

    
    def _MODEL(self, model):
        """ Machine learning model is created here. """
                    
        if self.algorithm == 'NN':
            _model = {'system': None,
                      'hiddenlayers': [6, 6],
                      'activation': ['Tanh', 'Tanh', 'Linear'],
                      'random_seed': None,
                      'epoch': 100,
                      'batch_size': None,
                      'device': 'cpu',
                      'force_coefficient': 0.03,
                      'stress_coefficient': None,
                      'stress_group': None,
                      'softmax_beta': None,
                      'alpha': None,
                      'unit': 'eV',
                      'restart': None,
                      'path': self.path,
                      'memory': 'in',
                      'optimizer': {},
                      }
            _model.update(model)
            
            if len(_model['activation']) != len(_model['hiddenlayers']) + 1:
                msg = '\nWarning: Incompatible activation functions and hiddenlayers.'
                print(msg)
                print('hiddenlayers: ', _model['hiddenlayers'])
                print('activations: ', _model['activation'])
                _model['activation'] = ['Tanh']*len(model['hiddenlayers'])+['Linear']
                print('revised activations: ', _model['activation'])
    
            if 'parameters' not in _model['optimizer']:
                _model['optimizer']['parameters'] = {}
            if 'derivative' not in _model['optimizer']:
                _model['optimizer']['derivative'] = True
            if 'method' not in _model['optimizer']:
                _model['optimizer']['method'] = 'lbfgs'

                # In full batch-LBFGS, epoch is 1. 
            if _model['optimizer']['method'] in ['lbfgs', 'LBFGS', 'lbfgsb']:
                if _model['batch_size'] is None: #full batch
                    _model['optimizer']['parameters']['max_iter'] = _model['epoch']
                    _model['epoch'] = 1
                else:
                    _model['optimizer']['parameters']['max_iter'] = 20

            self.model = NeuralNetwork(elements=_model['system'],
                                       hiddenlayers=_model['hiddenlayers'],
                                       activation=_model['activation'],
                                       random_seed=_model['random_seed'],
                                       epoch=_model['epoch'],
                                       batch_size=_model['batch_size'],
                                       device=_model['device'],
                                       alpha=_model['alpha'],
                                       force_coefficient=_model['force_coefficient'],
                                       stress_coefficient=_model['stress_coefficient'],
                                       stress_group=_model['stress_group'],
                                       softmax_beta=_model['softmax_beta'],
                                       unit=_model['unit'],
                                       restart=_model['restart'],
                                       path=_model['path'],
                                       memory=_model['memory'])
            self.optimizer = _model['optimizer']
                
        elif self.algorithm == 'PR':
            _model = {'system': None,
                      'force_coefficient': 0.0001,
                      'stress_coefficient': None,
                      'stress_group': None,
                      'order': 1,
                      'path': self.path,
                      'alpha': None,
                      'norm': 2,
                      'd_max': None,
                      }
            _model.update(model)
            self.model = PR(elements=_model['system'],
                            force_coefficient=_model['force_coefficient'],
                            stress_coefficient=_model['stress_coefficient'],
                            stress_group=_model['stress_group'],
                            order=_model['order'],
                            path=_model['path'],
                            alpha=_model['alpha'],
                            norm=_model['norm'],
                            d_max=_model['d_max'])
            self.optimizer = None


    def print_descriptors(self, _descriptors):
        """ Print the descriptors information. """

        print('Descriptor parameters:')
        keys = ['type', 'Rc', 'cutoff']
        for key in keys:
            print('{:12s}: {:}'.format(key, _descriptors[key]))

        if _descriptors['type'] in ['SO4', 'Bispectrum']:
            key_params = ['lmax', 'normalize_U']
        elif _descriptors['type'] in ['SO3', 'SOAP']:
            key_params = ['nmax', 'lmax', 'alpha']
        elif _descriptors['type'] in ['SNAP', 'snap']:
            key_params = ['lmax', 'rfac']
        elif _descriptors['type'] == 'EAD':
            key_params = ['L', 'eta', 'Rs']
        else:
            key_params = []

        for key in key_params:
            print('{:12s}: {:}'.format(key, _descriptors['parameters'][key]))
        print('\n')


    def print_logo(self):
        """ Print PyXtal_FF logo and version. """
        
        print("\n")
        print("""
               ______       _    _          _         _______ _______ 
              (_____ \     \ \  / /        | |       (_______|_______)
               _____) )   _ \ \/ / |_  ____| |        _____   _____   
              |  ____/ | | | )  (|  _)/ _  | |       |  ___) |  ___)  
              | |    | |_| |/ /\ \ |_( ( | | |_______| |     | |      
              |_|     \__  /_/  \_\___)_||_|_(_______)_|     |_|      
                     (____/      """)
        print("\n")
        print('        A Python package for Machine Learning Interatomic Force Field')
        print('         Developed by Zhu\'s group at University of Nevada Las Vegas')
        print('    The source code is available at https://github.com/qzhu2017/FF-project')
        print("\n")
        print('=========================== version', __version__,'=============================\n')


