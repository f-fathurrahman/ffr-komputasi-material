"""
Public API for XC ML models.
"""

from .base import XCBaseModel, ModelIOFiles, ModelKind, ModelSource, ConfigExt
from .calculator import MLXCCalculator

try:
    from .torch_backend import TorchXCModel  # noqa: F401
except ImportError:
    # Torch backend may not be available if torch is not installed
    pass

# Optional model definitions (require torch)
try:
    from .nn_models import ChannelEmbeddingMLP, ChannelEmbeddingResNet, MLP  # noqa: F401
except ImportError:
    # NN models may not be available if torch is not installed
    pass

__all__ = [
    "XCBaseModel",
    "ModelIOFiles",
    "ModelKind",
    "ModelSource",
    "ConfigExt",
    "MLXCCalculator",
    "TorchXCModel",
    "ChannelEmbeddingMLP",
    "ChannelEmbeddingResNet",
    "MLP",
]
