"""
Neural network model definitions for Delta Learning.

This module provides model architectures including:
- MLP: Simple Multi-Layer Perceptron
- ResidualBlock: Residual block with LayerNorm and SiLU activation
- ResidualBlockWithProjection: Residual block for dimension-changing layers
- ChannelEmbeddingResNet: ResNet-style architecture with channel embeddings
- ChannelEmbeddingMLP: MLP architecture with channel embeddings (no residual connections)
- FixedPBEEnhancementFactor: Fixed analytical PBE F_x(s) for debugging (no learnable params)
"""

from __future__ import annotations
from typing import List, Optional, TYPE_CHECKING, Dict, Any, Union

if TYPE_CHECKING:
    from ...data.data_manager import ComputationConfig


# Error messages
TORCH_NOT_AVAILABLE_ERROR = \
    "atom.xc.ml_xc.nn_models requires torch. " \
    "Install torch to use ChannelEmbeddingMLP and other NN models."

try:
    import torch
    import torch.nn as nn
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    torch = None
    nn = None


if not TORCH_AVAILABLE:
    raise ImportError(TORCH_NOT_AVAILABLE_ERROR)


if TORCH_AVAILABLE:
    
    class MLP(nn.Module):
        """Simple Multi-Layer Perceptron."""
        def __init__(self, input_dim=5, hidden_dims=[128, 64, 32], output_dim=1, dropout=0.0):
            super(MLP, self).__init__()
            layers = []
            prev_dim = input_dim
            
            for hidden_dim in hidden_dims:
                layers.append(nn.Linear(prev_dim, hidden_dim))
                layers.append(nn.Tanh())
                layers.append(nn.Dropout(dropout))
                prev_dim = hidden_dim
            
            layers.append(nn.Linear(prev_dim, output_dim))
            self.network = nn.Sequential(*layers)
        
        def forward(self, x):
            return self.network(x)
    
    
    class ResidualBlock(nn.Module):
        """Residual block with LayerNorm and SiLU activation."""
        def __init__(self, dim, dropout=0.0):
            super(ResidualBlock, self).__init__()
            self.ln = nn.LayerNorm(dim)
            self.block = nn.Sequential(
                nn.Linear(dim, dim),
                nn.SiLU(),
                nn.Dropout(dropout),
                nn.Linear(dim, dim),
                nn.SiLU(),
                nn.Dropout(dropout)
            )
        
        def forward(self, x):
            return x + self.block(self.ln(x))
    
    
    class ResidualBlockWithProjection(nn.Module):
        """Residual block for dimension-changing layers."""
        def __init__(self, in_dim, out_dim, dropout=0.0):
            super(ResidualBlockWithProjection, self).__init__()
            self.ln = nn.LayerNorm(in_dim)
            self.proj = nn.Linear(in_dim, out_dim)
            self.block = nn.Sequential(
                nn.Linear(in_dim, out_dim),
                nn.SiLU(),
                nn.Dropout(dropout),
                nn.Linear(out_dim, out_dim),
                nn.SiLU(),
                nn.Dropout(dropout)
            )
        
        def forward(self, x):
            return self.proj(x) + self.block(self.ln(x))
    
    
    class ChannelEmbeddingResNet(nn.Module):
        """
        Neural network with ResNet-style residual connections and LayerNorm.
        Each input channel goes through its own embedding network.
        """
        def __init__(self, input_dim=5, embedding_dims=[32, 16], hidden_dims=[128, 64, 32], 
                     output_dim=1, dropout=0.0, num_res_blocks=1):
            super(ChannelEmbeddingResNet, self).__init__()
            self.input_dim = input_dim
            
            # Create embedding networks for each channel
            self.channel_embeddings = nn.ModuleList()
            for i in range(input_dim):
                embedding_layers = []
                prev_dim = 1
                
                for emb_dim in embedding_dims:
                    embedding_layers.append(nn.Linear(prev_dim, emb_dim))
                    embedding_layers.append(nn.LayerNorm(emb_dim))
                    embedding_layers.append(nn.SiLU())
                    embedding_layers.append(nn.Dropout(dropout))
                    prev_dim = emb_dim
                
                final_embedding_dim = embedding_dims[-1] if embedding_dims else 1
                self.channel_embeddings.append(nn.Sequential(*embedding_layers))
            
            # Initial projection
            stacked_dim = input_dim * (embedding_dims[-1] if embedding_dims else 1)
            first_hidden_dim = hidden_dims[0] if hidden_dims else stacked_dim
            
            self.initial_proj = nn.Linear(stacked_dim, first_hidden_dim)
            self.initial_ln = nn.LayerNorm(first_hidden_dim)
            
            # Residual layers
            self.res_layers = nn.ModuleList()
            prev_dim = first_hidden_dim
            
            for hidden_dim in hidden_dims[1:]:
                if prev_dim != hidden_dim:
                    self.res_layers.append(ResidualBlockWithProjection(prev_dim, hidden_dim, dropout=dropout))
                else:
                    self.res_layers.append(ResidualBlock(prev_dim, dropout=dropout))
                prev_dim = hidden_dim
            
            for _ in range(num_res_blocks):
                if prev_dim > 0:
                    self.res_layers.append(ResidualBlock(prev_dim, dropout=dropout))
            
            final_dim = hidden_dims[-1] if hidden_dims else first_hidden_dim
            self.output_layer = nn.Linear(final_dim, output_dim)
        
        def forward(self, x):
            batch_size = x.size(0)
            
            embeddings = []
            for i in range(self.input_dim):
                channel = x[:, i:i+1]
                embedding = self.channel_embeddings[i](channel)
                embeddings.append(embedding)
            
            stacked = torch.cat(embeddings, dim=1)
            
            x = self.initial_proj(stacked)
            x = self.initial_ln(x)
            x = nn.SiLU()(x)
            
            for res_layer in self.res_layers:
                x = res_layer(x)
            
            output = self.output_layer(x)
            return output
    
    
    class ChannelEmbeddingMLP(nn.Module):
        """
        Neural network that embeds each input channel separately, then stacks and processes with MLP.
        Uses LayerNorm and SiLU activation, but NO residual connections (standard MLP).
        
        Architecture:
        1. Each input channel goes through its own embedding network (with LayerNorm)
        2. All embeddings are stacked together
        3. Stacked embeddings go through a standard MLP with LayerNorm and SiLU
        4. Final output layer
        """
        def __init__(self, input_dim=5, embedding_dims=[32, 16], hidden_dims=[128, 64, 32], 
                     output_dim=1, dropout=0.0):
            """
            Parameters
            ----------
            input_dim : int
                Number of input channels (features)
            embedding_dims : list
                Hidden dimensions for each channel's embedding network
            hidden_dims : list
                Hidden dimensions for the final MLP after stacking embeddings
            output_dim : int
                Output dimension
            dropout : float
                Dropout rate
            """
            super(ChannelEmbeddingMLP, self).__init__()
            self.input_dim = input_dim
            
            # Create embedding networks for each channel with LayerNorm
            self.channel_embeddings = nn.ModuleList()
            for i in range(input_dim):
                embedding_layers = []
                prev_dim = 1  # Each channel is a scalar
                
                for emb_dim in embedding_dims:
                    embedding_layers.append(nn.Linear(prev_dim, emb_dim))
                    embedding_layers.append(nn.LayerNorm(emb_dim))  # Add LayerNorm
                    embedding_layers.append(nn.SiLU())  # Replaced Tanh with SiLU
                    embedding_layers.append(nn.Dropout(dropout))
                    prev_dim = emb_dim
                
                # Final embedding dimension
                final_embedding_dim = embedding_dims[-1] if embedding_dims else 1
                self.channel_embeddings.append(nn.Sequential(*embedding_layers))
            
            # Final MLP after stacking all embeddings (standard MLP with LayerNorm and SiLU)
            stacked_dim = input_dim * (embedding_dims[-1] if embedding_dims else 1)
            mlp_layers = []
            prev_dim = stacked_dim
            
            for hidden_dim in hidden_dims:
                mlp_layers.append(nn.Linear(prev_dim, hidden_dim))
                mlp_layers.append(nn.LayerNorm(hidden_dim))  # Add LayerNorm
                mlp_layers.append(nn.SiLU())  # Replaced Tanh with SiLU
                mlp_layers.append(nn.Dropout(dropout))
                prev_dim = hidden_dim
            
            mlp_layers.append(nn.Linear(prev_dim, output_dim))
            self.final_mlp = nn.Sequential(*mlp_layers)
        
        def forward(self, x):
            """
            Forward pass.
            
            Parameters
            ----------
            x : torch.Tensor
                Input tensor of shape (batch_size, input_dim)
            
            Returns
            -------
            torch.Tensor
                Output tensor of shape (batch_size, output_dim)
            """
            # Split input into individual channels
            # x shape: (batch_size, input_dim)
            batch_size = x.size
            
            # Embed each channel separately
            embeddings = []
            for i in range(self.input_dim):
                # Extract i-th channel: (batch_size, 1)
                channel = x[:, i:i+1]
                # Pass through embedding network: (batch_size, embedding_dim)
                embedding = self.channel_embeddings[i](channel)
                embeddings.append(embedding)
            
            # Stack all embeddings: (batch_size, input_dim * embedding_dim)
            stacked = torch.cat(embeddings, dim=1)
            
            # Pass through final MLP (standard MLP, no residual)
            output = self.final_mlp(stacked)
            
            return output


    class FixedPBEEnhancementFactor(nn.Module):
        """
        Fixed analytical PBE exchange enhancement factor F_x(s) for debugging.
        
        Input: 2 channels [rho, grad_rho_reduced]. Only uses grad_rho_reduced (s).
        Output: scale * F_x(s) where F_x(s) = 1 + kappa - kappa/(1 + mu*s²/kappa).
        
        Has one learnable parameter `scale` (init=1). When target is PBE, training should
        converge with scale=1 and loss≈0. Use for pipeline debug.
        
        Usage:
            config = train_loader.get_computation_config()
            model = FixedPBEEnhancementFactor(
                features_list=["rho", "grad_rho_reduced"],
                config=config,
            )
            xc_model = TorchXCModel(
                model_kind="energy",
                features_list=["rho", "grad_rho_reduced"],
                model_cls=FixedPBEEnhancementFactor,
                model_init_kwargs={"features_list": ["rho", "grad_rho_reduced"], "config": config},
                ...
            )
        """
        EXPECTED_FEATURES = ["rho", "grad_rho_reduced"]
        KAPPA = 0.804
        MU = 0.2195149727645171  # PBE exchange

        def __init__(
            self,
            features_list: List[str],
            config: Union["ComputationConfig", Dict[str, Any]],
            model_kind: str = "energy",  # "energy" or "potential" - which output transform to use
            input_dim: int = 2,  # unused, for API compatibility
            output_dim: int = 1,  # unused, for API compatibility
            **kwargs,
        ):
            super().__init__()
            if features_list != self.EXPECTED_FEATURES:
                raise ValueError(
                    f"FixedPBEEnhancementFactor requires features_list={self.EXPECTED_FEATURES}, got {features_list}"
                )
            if model_kind not in ("energy", "potential"):
                raise ValueError(f"model_kind must be 'energy' or 'potential', got {model_kind}")
            self.features_list = list(features_list)
            self.config        = config
            self.model_kind    = model_kind
            self._s_idx        = features_list.index("grad_rho_reduced")
            # Learnable scale (init=1); when target is PBE, should converge to 1
            self.scale = nn.Parameter(torch.tensor(1.0, dtype=torch.float64))

        def _cfg_get(self, key: str, default=None):
            """Get config value; supports both dict and ComputationConfig."""
            cfg = self.config
            return cfg.get(key, default) if isinstance(cfg, dict) else getattr(cfg, key, default)

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            """
            x: features_transformed, shape (n_grid, 2)
            Returns: y_energy_transformed, shape (n_grid, 1)
            """
            from .torch_backend import TorchXCModel
            g = lambda k, d=None: self._cfg_get(k, d)
            # Inverse transform input: transformed -> physical
            features_physical = x.clone()
            if g("scale_features", True) and g("scaler_for_features") is not None:
                features_physical = TorchXCModel._torch_inverse_scaler(
                    features_physical, g("scaler_for_features")
                )
            if g("use_symlog_for_features", True) and g("linthresh_for_features") is not None:
                features_physical = TorchXCModel._torch_symexp(
                    features_physical, g("linthresh_for_features")
                )
            s = features_physical[:, self._s_idx]
            
            # PBE F_x(s) = 1 + kappa - kappa/(1 + mu*s²/kappa)
            kappa, mu = self.KAPPA, self.MU
            s2 = s * s  # avoid overflow in denominator
            F_x = 1.0 + kappa - kappa / (1.0 + (mu / kappa) * s2)
            out = self.scale * F_x.unsqueeze(-1)

            # Forward transform output: physical -> transformed (match target format)
            if self.model_kind == "potential":
                use_symlog = g("use_symlog_for_potential", True) and g("linthresh_for_potential") is not None
                use_scale  = g("scale_potential", True) and g("scaler_for_potential") is not None
                linthresh  = g("linthresh_for_potential")
                scaler     = g("scaler_for_potential")
            else:
                use_symlog = g("use_symlog_for_energy", True) and g("linthresh_for_energy") is not None
                use_scale  = g("scale_energy", True) and g("scaler_for_energy") is not None
                linthresh  = g("linthresh_for_energy")
                scaler     = g("scaler_for_energy")
                
            if use_symlog and linthresh is not None:
                out = TorchXCModel._torch_symlog(out, linthresh)
            if use_scale and scaler is not None:
                out = TorchXCModel._torch_apply_scaler(out, scaler)
            return out

