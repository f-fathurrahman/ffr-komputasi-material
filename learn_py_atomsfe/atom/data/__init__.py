"""Post-processing utilities for atomic DFT calculations.

This module provides:
- Data generation: DataGenerator, AtomicDataManager.generate_data
- Data loading: DataLoader, AtomicDataManager.load_data
- Data management: AtomicDataset, VxcDataLoader, AtomicDataManager
- Data processing: DataProcessor
- Model evaluation: evaluate_model, print_metrics (optional, requires sklearn)
- Visualization: DataVisualizer
"""


# Data generation utilities
from .data_generation import (  # noqa: F401
    TeeOutput,
    DataGenerator,
)

# Data loading utilities
from .data_loading import (  # noqa: F401
    DataLoader,
    SingleConfigurationData,
    FEATURE_ALIASES,
    NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL,
    NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY,
    VALID_FEATURES_LIST_FOR_POTENTIAL,
    VALID_FEATURES_LIST_FOR_ENERGY_DENSITY,
)

# Data manager
from .data_manager import (  # noqa: F401
    AtomicDataManager,
    AtomicDataset,
    VxcDataLoader,
)
DataManager = AtomicDataManager  # alias for convenience

# Data processing utilities
from .data_processing import (  # noqa: F401
    DataProcessor,
    VALID_SMOOTH_METHODS,
)

# Visualization utilities
from .data_visualization import (  # noqa: F401
    DataVisualizer,
)

__all__ = [
    # Data generation
    'TeeOutput',
    'DataGenerator',

    # Data loading
    'DataLoader',
    'SingleConfigurationData',
    'FEATURE_ALIASES',
    'NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL',
    'NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY',
    'VALID_FEATURES_LIST_FOR_POTENTIAL',
    'VALID_FEATURES_LIST_FOR_ENERGY_DENSITY',

    # Data manager
    'AtomicDataManager',
    'AtomicDataset',
    'VxcDataLoader',
    'DataManager',

    # Data processing
    'DataProcessor',
    'VALID_SMOOTH_METHODS',

    # Visualization
    'DataVisualizer',
]
