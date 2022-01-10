from collections import Counter
import os
import warnings

import torch
from torch.nn import Module

from regression_gppotential import context_setting, AutoMean, _regression
from util_parallel import if_master
from GaussianProcessPotential import GaussianProcessPotential
from similarity_similarity import SimilarityKernel
from descriptor_atoms import Local, TorchAtoms, AtomsData, LocalsData
from util_util import mkdir_p
from regression_scores import coeff_of_determination


class PosteriorPotential(Module):

    def __init__(self, gp, data=None, inducing=None, group=None, **setting):
        #
        super().__init__()
        print("Initializing PosteriorPotential")
        #
        if type(gp) == GaussianProcessPotential:
            self.gp = gp
        elif type(gp) == list:
            self.gp = GaussianProcessPotential(gp)
        elif issubclass(gp.__class__, SimilarityKernel):
            self.gp = GaussianProcessPotential([gp])
        else:
            raise RuntimeError(f'type {type(gp)} is not recognized')
        #
        if group is not None:
            self.attach_process_group(group)
        else:
            self.is_distributed = False
        #
        self.ignore_forces = False
        #
        if data is not None:
            self.set_data(data, inducing=inducing, **setting)
        else:
            self.data = AtomsData([])
            self.X = LocalsData([])
            self.has_target_forces = False
            self.M = torch.empty(0, 0)
            self.Ke = torch.empty(0, 0)
            self.Kf = torch.empty(0, 0)

    @property
    def ndata(self):
        return len(self.data)

    @context_setting
    def set_data(self, _data, inducing=None):
        data = _data.subset(self.gp.species)
        self.data = data
        if inducing is None:
            raise RuntimeWarning(
                'This (inducing=None) has not been used in long while!')
            p = self.gp(data, inducing)
            self.X = copy.deepcopy(data)  # TODO: consider not copying
            self.mu = p.precision_matrix @ (gp.Y(data)-p.loc)
            self.nu = p.precision_matrix
            self.has_target_forces = True
        else:
            X = inducing.subset(self.gp.species)
            self.Ke = self.gp.kern(data, X, cov='energy_energy')
            self.Kf = self.gp.kern(data, X, cov='forces_energy')
            if data.is_distributed:
                torch.distributed.all_reduce(self.Ke)
                torch.distributed.all_reduce(self.Kf)
            self.M = self.gp.kern(X, X, cov='energy_energy')
            self.X = X
            self.make_munu()
            self.has_target_forces = False

    @property
    def noise(self):
        return self.gp.noise.signal

    @noise.setter
    def noise(self, value):
        self.gp.noise.signal = value
        self.make_munu()

    @property
    def descriptors(self):
        return self.gp.descriptors

    @property
    def cutoff(self):
        return self.gp.cutoff

    @property
    def _cutoff(self):
        return self.cutoff

    @property
    def inducing(self):
        return self.X

    @property
    def K(self):
        return torch.cat([self.Ke, self.Kf], dim=0)

    @property
    def mean(self):
        return self.gp.parametric

    @mean.setter
    def mean(self, value):
        self.gp.parametric = value

    def make_munu(self, algo=2, noisegrad=False, **kw):
        if self.M.numel() == 0 or self.K.numel() == 0:
            return
        if not isinstance(self.mean, AutoMean):
            self.mean = AutoMean()  # xxxMean -> AutoMean
        self.mean.set_data(self.data)
        parallel = torch.distributed.is_initialized()
        if parallel:
            rank = torch.distributed.get_rank()
        else:
            rank = 0
        if rank == 0:
            if algo == 0:
                # allocates too much memory
                self.mu, self.nu, self.ridge, self.choli = projected_process_auxiliary_matrices_D(
                    self.K, self.M, self.gp.Y(self.data), self.gp.diagonal_ridge(self.data), chol_inverse=True)
            elif algo == 1:
                L, ridge = jitcholesky(self.M)
                self.ridge = torch.as_tensor(ridge)
                #sigma = self.gp.diagonal_ridge(self.data).sqrt()
                sigma = self.gp.noise.signal*self.M.diag().mean()
                if not noisegrad:
                    sigma = sigma.detach()
                #A = torch.cat((self.K/sigma.view(-1, 1), L.t()))
                A = torch.cat((self.K, sigma.view(1)*L.t()))
                #Y = torch.cat((self.gp.Y(self.data)/sigma, torch.zeros(L.size(0))))
                Y = torch.cat((self.gp.Y(self.data), torch.zeros(L.size(0))))
                Q, R = torch.linalg.qr(A)
                self.mu = (R.inverse()@Q.t()@Y).contiguous()
                # self.nu = None # is not needed anymore
                self.choli = L.inverse().contiguous()
            elif algo == 2:
                _regression(self, optimize=False)
            elif algo == 3:
                _regression(self, optimize=True, **kw)
        else:
            self.ridge = torch.zeros([])
            self.mu = torch.zeros_like(self.M[0])
            self.choli = torch.zeros_like(self.M)
        if parallel:
            torch.distributed.broadcast(self.ridge, 0)
            torch.distributed.broadcast(self.mu, 0)
            torch.distributed.broadcast(self.choli, 0)
            self.mean.sync_params()
        if not noisegrad and (self.mu.requires_grad or self.choli.requires_grad):
            warnings.warn('mu or choli requires grad!')
        self.Mi = self.choli.t()@self.choli
        self.make_stats()

    def optimize_model_parameters(self, **kw):
        self.make_munu(algo=3, **kw)

    def make_stats(self):
        n = len(self.data)
        y = self.gp.Y(self.data)
        yy = self.K@self.mu
        diff = yy - y
        self._ediff = diff[:n]/torch.tensor(self.data.natoms)
        self._fdiff = diff[n:]
        self._force_r2 = coeff_of_determination(yy[n:], y[n:])
        self._stats = [self._ediff.mean(), self._ediff.abs().mean(),
                       self._fdiff.mean(), self._fdiff.abs().mean(),
                       self._force_r2]
        # forces info
        self._f_max = y[n:].abs().max()
        self._f_std = y[n:].var().sqrt()
        # needed for special cases
        self.indu_counts = Counter()
        self.kern_diag_mean = Counter()
        for i, x in enumerate(self.X):
            self.indu_counts[x.number] += 1
            self.kern_diag_mean[x.number] += float(self.M[i, i])
        for num in self.kern_diag_mean.keys():
            self.kern_diag_mean[num] /= self.indu_counts[num]
        #
        # predictive variance for x
        # _vscale *[ k(x,x) - k(x,m)k(m,m)^{-1}k(m,x) ]
        self.indu_numbers = torch.tensor([x.number for x in self.X])
        self._vscale = {}
        mu = self.mu*(self.M@self.mu)
        for z in self.indu_counts.keys():
            I = self.indu_numbers == z
            self._vscale[z] = mu[I].sum()/I.sum()

    @property
    def sigma_e(self):
        return self._stats[1]

    @property
    def sigma_f(self):
        return self._stats[3]

    def is_ok(self):
        e_ok = (self._stats[0]-self._stats[1]) * \
            (self._stats[0]+self._stats[1]) < 0
        f_ok = (self._stats[2]-self._stats[3]) * \
            (self._stats[2]+self._stats[3]) < 0
        return e_ok and f_ok

    def is_well(self, a=None, b=None):
        x = True if a is None else abs(self._stats[0]) < a*self._stats[1]
        y = True if b is None else abs(self._stats[2]) < b*self._stats[3]
        return all([self.is_ok(), x, y])

    def tune_noise(self, a=None, b=None, lr=1., min_steps=0, weighted=lambda a: a, verbose=False):

        def step():
            opt.zero_grad()
            self.make_munu(noisegrad=True)
            losses = -torch.distributions.normal.Normal(
                0., self._stats[1]).log_prob(self._ediff)
            loss = weighted(losses).sum()
            loss.backward()
            opt.step()
            return loss

        noise = list(self.gp.noise.parameters())[0]
        noise.requires_grad = True
        opt = torch.optim.Adam([noise], lr=lr)
        steps = 0
        _min = None
        while not self.is_well(a, b) or steps < min_steps:
            loss = step()
            steps += 1
            if (_min is None) or loss < _min:
                _min = loss
                _min_arg = self.gp.noise.signal
                _min_step = steps
            if verbose:
                print(
                    f'{steps}  loss: {loss}  noise: {self.gp.noise.signal}  global: ({_min_step})')
        self.make_munu(noisegrad=False)
        return steps, _min, _min_arg

    @property
    def ref_M(self):
        return self.M + self.ridge*torch.eye(self.M.size(0))

    @context_setting
    def leakage(self, loc):
        a = self.gp.kern(self.X, loc, cov='energy_energy')
        b = self.choli @ a
        c = b.t()@b
        d = self.gp.kern(loc, loc, cov='energy_energy') + self.ridge
        return (1-c/d).view(1)

    def leakages(self, X):
        return torch.cat([self.leakage(x) for x in iterable(X)])

    @context_setting
    def remake_all(self):
        self.set_data(self.data, self.X)

    def add_kernels(self, kernels, remake_all=True):
        self.gp.add_kernels(kernels)
        self.data.apply('add_descriptors', kernels, dont_save_grads=False)
        self.inducing.stage(kernels, dont_save_grads=True)
        if remake_all:
            self.remake_all()

    @context_setting
    def add_data(self, data, remake=True):
        assert data[0].includes_species(self.gp.species)
        Ke = self.gp.kern(data, self.X, cov='energy_energy')
        Kf = self.gp.kern(data, self.X, cov='forces_energy')
        if (data[0].is_distributed if type(data) == list else data.is_distributed):
            torch.distributed.all_reduce(Ke)
            torch.distributed.all_reduce(Kf)
        self.Ke = torch.cat([self.Ke, Ke], dim=0)
        self.Kf = torch.cat([self.Kf, Kf], dim=0)
        self.data += data
        if remake:
            self.make_munu()

    @context_setting
    def add_inducing(self, X, col=None, remake=True):
        assert X.number in self.gp.species
        Ke = self.gp.kern(self.data, X, cov='energy_energy')
        Kf = self.gp.kern(self.data, X, cov='forces_energy')
        if self.data.is_distributed:
            torch.distributed.all_reduce(Ke)
            torch.distributed.all_reduce(Kf)
        if self.Ke.numel() > 0:
            self.Ke = torch.cat([self.Ke, Ke], dim=1)
            self.Kf = torch.cat([self.Kf, Kf], dim=1)
        else:
            self.Ke = Ke
            self.Kf = Kf
        if col is None:
            a = self.gp.kern(self.X, X, cov='energy_energy')
        else:
            a = col
        b = self.gp.kern(X, X, cov='energy_energy')
        self.M = torch.cat(
            [torch.cat([self.M, a.t()]), torch.cat([a, b])], dim=1)
        self.X += X
        if remake:
            self.make_munu()

    def pop_1data(self, remake=True, clear_cached=True):
        self.Ke = self.Ke[:-1]
        self.Kf = self.Kf[:-3*self.data[-1].natoms]
        if clear_cached:
            self.gp.clear_cached([self.data.X[-1]])
        del self.data.X[-1]
        if remake:
            self.make_munu()

    def pop_1inducing(self, remake=True, clear_cached=True):
        self.Ke = self.Ke[:, :-1]
        self.Kf = self.Kf[:, :-1]
        self.M = self.M[:-1, :-1]
        if clear_cached:
            self.gp.clear_cached([self.X.X[-1]])
        del self.X.X[-1]
        if remake:
            self.make_munu()

    def popfirst_1data(self, remake=True, clear_cached=True):
        self.Ke = self.Ke[1:]
        self.Kf = self.Kf[3*self.data[0].natoms:]
        if clear_cached:
            self.gp.clear_cached([self.data.X[0]])
        del self.data.X[0]
        if remake:
            self.make_munu()

    def popfirst_1inducing(self, remake=True, clear_cached=True):
        self.Ke = self.Ke[:, 1:]
        self.Kf = self.Kf[:, 1:]
        self.M = self.M[1:, 1:]
        if clear_cached:
            self.gp.clear_cached([self.X.X[0]])
        del self.X.X[0]
        if remake:
            self.make_munu()

    def downsize(self, n, m, first=False, lii=False, remake=True):
        """
        lii: remove least important inducing
        """
        if any([x.requires_grad for x in (self.M, self.Ke, self.Kf)]):
            raise RuntimeError('cov matrices require grad!')
        ch1 = 0
        while len(self.data) > n:
            if first:
                self.popfirst_1data(remake=False)
            else:
                self.pop_1data(remake=False)
            ch1 += 1
        ch2 = 0
        if lii and m < len(self.X):
            indices = torch.argsort(self.M.sum(axis=1)).tolist()
            ch2 = indices[:m]
            self.select_inducing(ch2, deleted=indices[m:], remake=False)
        else:
            while len(self.X) > m:
                if first:
                    self.popfirst_1inducing(remake=False)
                else:
                    self.pop_1inducing(remake=False)
                ch2 += 1
        if remake and (ch1 or ch2):
            self.make_munu()
        return ch1, ch2

    def add_1atoms(self, atoms, ediff, fdiff):
        """ ediff here is ediff_tot """
        if not atoms.includes_species(self.gp.species):
            return 0, 0, 0
        kwargs = {'use_caching': True}
        #
        if len(self.data) == 0:
            if len(self.X) > 0:
                self.add_data([atoms], **kwargs)
            else:
                self.data.append(atoms)
            return 1, float('inf'), float('inf')
        #
        use_forces = fdiff < float('inf') and not self.ignore_forces
        #
        e1 = self([atoms], all_reduce=atoms.is_distributed, **kwargs)
        if use_forces:
            f1 = self([atoms], 'forces',
                      all_reduce=atoms.is_distributed, **kwargs)
        self.add_data([atoms], **kwargs)
        e2 = self([atoms], all_reduce=atoms.is_distributed, **kwargs)
        if use_forces:
            f2 = self([atoms], 'forces',
                      all_reduce=atoms.is_distributed, **kwargs)
        #
        de = abs(e1-e2)
        df = 0.
        if not use_forces:
            reject = de < ediff
        else:
            # TODO: better algorithm!
            d = (f2-f1).view(-1)
            df = d.abs().mean()
            df_max = d.abs().max()
            #reject = de < ediff and df < fdiff and df_max < 3*fdiff
            N = torch.distributions.Normal(0., fdiff)
            reject = de < ediff and N.log_prob(d).mean() > N.log_prob(fdiff)
        #
        blind = torch.cat([e1, e2]).allclose(torch.zeros(1))
        if reject and not blind:
            self.pop_1data(clear_cached=True)
            added = 0
        else:
            added = 1
        return added, de, df

    def add_1atoms_fast(self, atoms, ediff, fdiff, xyz, cov, is_distributed):
        """ ediff here is ediff_tot """
        if not atoms.includes_species(self.gp.species):
            return 0, 0, 0
        kwargs = {'use_caching': True}
        #
        if len(self.data) == 0:
            if len(self.X) > 0:
                self.add_data([atoms], **kwargs)
            else:
                self.data.append(atoms)
            return 1, float('inf'), float('inf')
        #
        use_forces = fdiff < float('inf') and not self.ignore_forces
        f1 = torch.zeros_like(xyz)
        f2 = torch.zeros_like(xyz)
        #
        e1 = (cov@self.mu).sum().view(1)
        if e1.grad_fn and use_forces:
            f1 = -torch.autograd.grad(e1, xyz, retain_graph=True)[0]
        self.add_data([atoms], **kwargs)
        e2 = (cov@self.mu).sum().view(1)
        if e2.grad_fn and use_forces:
            f2 = -torch.autograd.grad(e2, xyz, retain_graph=True)[0]
        if is_distributed:
            torch.distributed.all_reduce(e1)
            torch.distributed.all_reduce(e2)
            if use_forces:
                torch.distributed.all_reduce(f1)
                torch.distributed.all_reduce(f2)
        #
        de = abs(e1-e2)
        df = 0.
        if not use_forces:
            reject = de < ediff
        else:
            # TODO: better algorithm!
            d = (f2-f1).view(-1)
            df = d.abs().mean()
            df_max = d.abs().max()
            fdiff_t = torch.tensor(fdiff)
            #reject = de < ediff and df < fdiff and df_max < 3*fdiff
            N = torch.distributions.Normal(0., fdiff_t)
            # reject = de < ediff and N.log_prob(d).mean() > N.log_prob(fdiff_t)
            reject = (N.log_prob(d).mean() > N.log_prob(fdiff_t) and
                      df_max < 3*fdiff)
        #
        blind = torch.cat([e1, e2]).allclose(torch.zeros(1))
        if reject and not blind:
            self.pop_1data(clear_cached=True)
            added = 0
        else:
            added = 1
        return added, de, df

    def add_1inducing(self, _loc, ediff, detach=True):
        if _loc.number not in self.gp.species:
            return 0, 0.
        kwargs = {'use_caching': True}
        if detach:
            loc = _loc.detach()
            loc.stage(self.gp.kern.kernels, dont_save_grads=True)
        else:
            loc = _loc
        #
        if len(self.X) == 0:
            if len(self.data) > 0:
                self.add_inducing(loc, **kwargs)
            else:
                self.X += loc
            return 1, float('inf')
        #
        e1 = self(loc, **kwargs)
        self.add_inducing(loc, **kwargs)
        e2 = self(loc, **kwargs)
        de = abs(e1-e2)
        blind = torch.cat([e1, e2]).allclose(torch.zeros(1))
        if (de < ediff and not blind) or self.ridge > 0.:
            self.pop_1inducing(clear_cached=True)
            added = 0
        else:
            added = 1
        return added, de

    def add_ninducing(self, _locs, ediff, detach=True, descending=True, leaks=None):
        selected = torch.as_tensor([i for i, loc in enumerate(_locs)
                                    if loc.number in self.gp.species])
        locs = [_locs[i] for i in selected]
        if len(locs) == 0:
            return 0, 0.
        if descending:
            if leaks is None:
                _leaks = self.leakages(locs)
            else:
                _leaks = leaks.index_select(0, selected)
            q = torch.argsort(_leaks, descending=True)
        else:
            q = torch.arange(len(locs))
        added_refs = 0
        for k in q:
            loc = locs[k]
            _ediff = ediff if len(self.X) > 1 else torch.finfo().eps
            added, change = self.add_1inducing(loc, _ediff, detach=detach)
            if added:
                added_refs += 1
            else:
                if descending:
                    break
        return added_refs, change

    def as_(self, _obj, group=None):
        if type(_obj) == Local:
            obj = _obj
            obj.stage(self.descriptors)
        elif type(_obj) == TorchAtoms:
            obj = _obj
        else:
            if group is None and torch.distributed.is_initialized():
                group = torch.distributed.group.WORLD
            obj = TorchAtoms(ase_atoms=_obj, cutoff=self.cutoff,
                             descriptors=self.descriptors, group=group)
        return obj

    def eat(self, _atoms, ediff, fdiff, group=None):
        atoms = self.as_(_atoms, group=group)
        if len(self.data) == 0:
            i = atoms.first_of_each_atom_type()
            locs = atoms.gathered()
            inducing = LocalsData([locs[j] for j in i])
            data = AtomsData([atoms])
            self.set_data(data, inducing)
            remaining = [locs[j] for j in range(atoms.natoms) if j not in i]
            self.add_ninducing(remaining, ediff)
        else:
            if atoms.is_distributed:
                leaks = torch.zeros(atoms.natoms)
                leaks[atoms.indices] = self.leakages(atoms.loc)
                torch.distributed.all_reduce(leaks)
            else:
                leaks = None
            locs = atoms.gathered()
            added_refs, change = self.add_ninducing(locs, ediff, leaks=leaks)
            if added_refs > 0:
                self.add_1atoms(atoms, ediff, fdiff)

    def select_inducing(self, indices, deleted=None, remake=True):
        i = torch.as_tensor(indices)
        self.Ke = self.Ke.index_select(1, i)
        self.Kf = self.Kf.index_select(1, i)
        self.M = self.M.index_select(0, i).index_select(1, i)
        if deleted:
            self.gp.clear_cached([self.X.X[j] for j in deleted])
        self.X.X = [self.X.X[j] for j in i]
        if remake:
            self.make_munu()

    def attach_process_group(self, *args, **kwargs):
        self.gp.attach_process_group(*args, **kwargs)
        self.is_distributed = True

    def detach_process_group(self, *args, **kwargs):
        self.gp.detach_process_group(*args, **kwargs)
        self.is_distributed = False

    def train(self, *args, **kwargs):
        train_gpp(self.gp, *args, **kwargs)

    def save(self, file, supress_warnings=True):
        cached = self.gp.cached
        self.gp.del_cached()
        data = self.data
        del self.data
        self._raw_data = [atoms.as_ase() for atoms in data]
        with warnings.catch_warnings():
            if supress_warnings:
                warnings.simplefilter("ignore")
            torch.save(self, file)
        self.data = data
        self.gp.cahced = cached

    @if_master
    def to_folder(self, folder, info=None,
                  overwrite=True, supress_warnings=True, pickle_data=False,
                  to_traj=False):

        if pickle_data and self.data.is_distributed:
            raise NotImplementedError(
                'trying to pickle data which is distributed! call gathere_() first!')
        if not overwrite:
            folder = safe_dirname(folder)
        
        mkdir_p(folder)
        
        with open(os.path.join(folder, 'cutoff'), 'w') as file:
            file.write('{}\n'.format(self.cutoff))

        if to_traj:  # not necessary, data will be pickled as self._raw_data
            self.data.to_traj(os.path.join(folder, 'data.traj'))
            self.X.to_traj(os.path.join(folder, 'inducing.traj'))

        self.gp.to_file(os.path.join(folder, 'gp'))

        self.save( os.path.join(folder, 'model'), supress_warnings=supress_warnings )

        # pickles (inducing are pickled with model)
        if pickle_data:
            with warnings.catch_warnings():
                if supress_warnings:
                    warnings.simplefilter("ignore")
                torch.save(self.data, os.path.join(folder, 'data.pckl'))

        # info
        with open(os.path.join(folder, 'info'), 'w') as file:
            file.write('data: {}, inducing: {}\n'.format(
                len(self.data), len(self.X)))
            if info is not None:
                if type(info) == str:
                    file.write('{}\n'.format(info))
                elif hasattr(info, '__iter__'):
                    for inf in info:
                        file.write('{}\n'.format(inf))
                else:
                    file.write('{}\n'.format(info))

        # stats
        with open(os.path.join(folder, 'stats'), 'w') as file:
            e1, e2, f1, f2, r2 = (float(v) for v in self._stats)
            file.write(f'ediff -> mean: {e1} std: {e2}  ')
            file.write(f'fdiff -> mean: {f1} std: {f2}  ')
            file.write(f'R2: {r2}\n')

    @context_setting
    def forward(self, test, quant='energy', variance=False, all_reduce=False):
        #
        print("Pass here in PosteriorPotential 1043")
        #
        shape = {'energy': (-1,), 'forces': (-1, 3)}
        A = self.gp.kern(test, self.X, cov=quant+'_energy')
        if self.has_target_forces:
            A = torch.cat([A, self.gp.kern(test, self.X, cov=quant+'_forces')],
                          dim=1)
        out = (A @ self.mu).view(*shape[quant])
        if all_reduce:
            torch.distributed.all_reduce(out)
        # mean is not distributed!
        if quant == 'energy':
            mean = self.mean(test, forces=False)
        else:
            _, mean = self.mean(test, forces=True)
        out = out + mean
        if variance:
            if all_reduce:
                raise NotImplementedError(
                    'all_reduce with variance=True is not implemented')
            var = (self.gp.kern.diag(test, quant) -
                   (A @ self.nu @ A.t()).diag()).view(*shape[quant])
            return out, var
        else:
            return out

    @context_setting
    def predict(self, test, variance=False):
        if self.gp.parametric is not None:
            raise NotImplementedError(
                'this method is not updated to include parametric potential')
        A = self.gp.kern(test, self.X, cov='energy_energy')
        B = self.gp.kern(test, self.X, cov='forces_energy')
        if self.has_target_forces:
            A = torch.cat([A, self.gp.kern(test, self.X, cov='energy_forces')],
                          dim=1)
            B = torch.cat([B, self.gp.kern(test, self.X, cov='forces_forces')],
                          dim=1)

        em, fm = self.mean(test, forces=True)
        energy = A @ self.mu + em
        forces = B @ self.mu + fm

        out = (energy, forces.view(-1, 3))

        if variance:
            energy_var = (self.gp.kern.diag(test, 'energy') -
                          (A @ self.nu @ A.t()).diag())
            forces_var = (self.gp.kern.diag(test, 'forces') -
                          (B @ self.nu @ B.t()).diag())
            out += (energy_var, forces_var.view(-1, 3))

        return out

def PosteriorPotentialFromFolder(folder, load_data=True, update_data=True, group=None, distrib=None):
    from descriptor_atoms import AtomsData
    from util_caching import strip_uid
    self = torch.load(os.path.join(folder, 'model'))
    strip_uid(self.X)
    if load_data:
        if os.path.isfile(os.path.join(folder, 'data.pckl')):
            self.data = torch.load(os.path.join(folder, 'data.pckl'))
            strip_uid(self.data)
            if group:
                self.data.distribute_(group)
        else:
            if hasattr(self, '_raw_data'):
                self.data = AtomsData(
                    self._raw_data, convert=True, group=group, ranks=distrib)
                del self._raw_data
            else:  # for backward compatibility
                self.data = AtomsData(traj=os.path.join(folder, 'data.traj'),
                                      group=group, ranks=distrib)
            if update_data:
                self.data.update(
                    cutoff=self.cutoff, descriptors=self.gp.kern.kernels)
    return self

