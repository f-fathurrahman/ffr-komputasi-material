import math
from typing import Sequence, Tuple

import numpy as np

import torch
from torch import Tensor
from torch import nn
import pytorch_lightning as pl


def load_data():
    coord = np.array(np.load("../DATASET/DeepPot/input_coord.npy"), dtype="float32")
    atom_types = np.loadtxt("../DATASET/DeepPot/type.raw", dtype=int)
    elems = np.unique(atom_types).tolist()
    atom_types = np.array([[elems.index(i) for i in atom_types]])
    atom_types = atom_types.repeat(len(coord), axis=0)
    grad = np.array(np.load("../DATASET/DeepPot/input_grad.npy"), dtype="float32")

    # Now, convert them to Tensor
    device = torch.device("cpu")
    coord = torch.from_numpy(coord).to(device)
    atom_types = torch.from_numpy(atom_types).to(device)
    grad = torch.from_numpy(grad).to(device)
    return coord, atom_types, grad



class MySequential(nn.Sequential):
    def forward(self, input: Tuple[Tensor, Tensor]) -> Tuple[Tensor, Tensor]:
        for module in self:
            input = module(input)
        return input

class MyDense(nn.Module):
    def __init__(self, num_channels: int, in_features: int, out_features: int, bias: bool = True, activation: bool = False, residual: bool = False) -> None:
        super().__init__()
        self.num_channels = num_channels
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(torch.Tensor(num_channels, out_features, in_features))
        if bias:
            self.bias = nn.Parameter(torch.Tensor(num_channels, out_features))
        else:
            self.register_parameter('bias', None)
        self.activation = activation
        self.residual = residual
        self.reset_parameters()

    def reset_parameters(self) -> None:
        for w in self.weight:
            nn.init.kaiming_uniform_(w, a=math.sqrt(5))
        if self.bias is not None:
            for b, w in zip(self.bias, self.weight):
                fan_in, _ = nn.init._calculate_fan_in_and_fan_out(w)
                bound = 1 / math.sqrt(fan_in)
                nn.init.uniform_(b, -bound, bound)

    def forward(self, input: Tuple[Tensor, Tensor]) -> Tuple[Tensor, Tensor]:
        x, channels = input
        weight: Tensor = self.weight[channels]
        output: Tensor = torch.bmm(x.transpose(0, 1), weight.transpose(1, 2)).transpose(0, 1)

        if self.bias is not None:
            bias = self.bias[channels]
            output = output + bias

        if self.activation:
            output = torch.tanh(output)

        if self.residual:
            if output.shape[2] == x.shape[2]:
                output = output + x
            elif output.shape[2] == x.shape[2] * 2:
                output = output + torch.cat([x, x], dim=2)
            else:
                raise NotImplementedError("Not implemented")

        return output, channels

    def extra_repr(self) -> str:
        return 'num_channels={}, in_features={}, out_features={}, bias={}, activation={}, residual={}'.format(
            self.num_channels, self.in_features, self.out_features, self.bias is not None, self.activation, self.residual
        )



def local_environment(coords: Tensor) -> Tuple[Tensor, Tensor]:
    num_batches, num_channels, _ = coords.size()
    rij = coords[:, :, None] - coords[:, None]
    dij = torch.norm(rij, dim=3)
    mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device) # remove self-interaction
    rij = torch.masked_select(rij, mask.unsqueeze(2)).view(num_batches, num_channels, num_channels - 1, 3)
    dij = torch.masked_select(dij, mask).view(num_batches, num_channels, num_channels - 1)
    dij_inv = 1 / dij
    dij2_inv = dij_inv * dij_inv

    loc_env_r = dij_inv
    loc_env_a = rij * dij2_inv.unsqueeze(3)

    return loc_env_r, loc_env_a



class Feature(nn.Module):
    def __init__(self, n_types: int, neuron: Sequence[int] = [25, 50, 100], axis_neuron: int = 4) -> None:
        super().__init__()
        self.n_types = n_types
        self.neuron = neuron
        self.axis_neuron = axis_neuron

        layers = [MyDense(n_types * n_types, 1, neuron[0], activation=True)]
        for i in range(len(neuron)-1):
            layers.append(MyDense(n_types * n_types, neuron[i], neuron[i+1], activation=True, residual=True))
        self.local_embedding = MySequential(*layers)

    def forward(self, coords: Tensor, atom_types: Tensor) -> Tensor:
        num_batches, num_channels, _ = coords.size()
        loc_env_r, loc_env_a = local_environment(coords)

        neighbor_types = atom_types.repeat(num_channels, 1)
        mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device)
        neighbor_types = torch.masked_select(neighbor_types, mask).view(num_channels, -1)
        indices = ((atom_types * self.n_types).unsqueeze(-1) + neighbor_types).view(-1)

        output, _ = self.local_embedding((loc_env_r.view(num_batches, -1, 1), indices))
        output = output.view(num_batches, num_channels, num_channels - 1, -1)

        output = torch.transpose(output, 2, 3) @ (loc_env_a @ (torch.transpose(loc_env_a, 2, 3) @ output[..., :self.axis_neuron]))
        output = output.view(num_batches, num_channels, -1)

        return output

    @property
    def output_length(self) -> int:
        return self.neuron[-1] * self.axis_neuron


class Fitting(nn.Module):
    def __init__(self, n_types: int, in_features: int, neuron: Sequence[int] = [240, 240, 240]) -> None:
        super().__init__()
        layers = [MyDense(n_types, in_features, neuron[0], activation=True)]
        for i in range(len(neuron)-1):
            layers.append(MyDense(n_types, neuron[i], neuron[i+1], activation=True, residual=True))
        layers.append(MyDense(n_types, neuron[-1], 1))
        self.fitting_net = MySequential(*layers)

    def forward(self, input : Tuple[Tensor, Tensor]) -> Tensor:
        output, _ = self.fitting_net(input)
        return output



class DeepPot(pl.LightningModule):
    def __init__(self, descriptor: nn.Module, fitting_net: nn.Module, learning_rate=5e-4) -> None:
        super().__init__()
        self.descriptor = descriptor
        self.fitting_net = fitting_net
        self.learning_rate = learning_rate

    def forward(self, coords: torch.Tensor, atom_types: torch.Tensor):
        coords.requires_grad_()
        descriptors = self.descriptor(coords, atom_types) # call descriptor
        atomic_energies = self.fitting_net((descriptors, atom_types)) # evaluate
        energy = torch.unbind(torch.sum(atomic_energies, dim=1))
        gradient, = torch.autograd.grad(energy, [coords], create_graph=True)
        return torch.hstack(energy), gradient

    def training_step(self, batch, batch_idx):
        qm_coord, atom_types, grad = batch
        ene_pred, grad_pred = self(qm_coord, atom_types[0])
        loss = nn.functional.mse_loss(grad_pred, grad)
        self.log('train_loss', loss)
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        scheduler = {'scheduler': torch.optim.lr_scheduler.ExponentialLR(optimizer, 0.95),
                     'interval': 'epoch',
                     'frequency': 1,
                    }
        return [optimizer], [scheduler]


# Model
descriptor = Feature(4, neuron=[25, 50], axis_neuron=4)
fitting_net = Fitting(4, descriptor.output_length, neuron=[120, 120])
model = DeepPot(descriptor, fitting_net, learning_rate=5e-4)

# Data
coord, atom_types, grad = load_data()

# forward pass
res = model(coord[0:1,:,:], atom_types[0]) # xxx only need one instance of atom_types?


