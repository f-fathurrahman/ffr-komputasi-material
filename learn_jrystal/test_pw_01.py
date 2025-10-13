import jax
import jax.numpy as jnp
import my_jrystal as jr
import numpy as np

jax.config.update("jax_enable_x64", True)

key = jax.random.PRNGKey(123)
diamond_file_path = "./geometry/diamond.xyz"
crystal = jr.Crystal.create_from_file(diamond_file_path)
num_bands = crystal.num_electron
key = jax.random.PRNGKey(123)
kpts = jr.grid.k_vectors(crystal.A, [1, 1, 1])
g_vecs = jr.grid.g_vectors(crystal.A, [7, 8, 9])
freq_mask = jr.grid.cubic_mask([7, 8, 9])


def test_wave():
    pw_param = jr.pw.param_init(key, num_bands, 1, freq_mask)
    coeff = jr.pw.coeff(pw_param, freq_mask)
    wave_grid1 = jr.pw.wave_grid(coeff, crystal.vol)

    @jr.utils.vmapstack(3)
    def wave(r):
        return jr.pw.wave_r(r, coeff, crystal.A, g_vecs)

    r_vecs = jr.grid.g2r_vector_grid(g_vecs, crystal.A)
    wave_grid2 = wave(r_vecs)
    wave_grid2 = jnp.transpose(wave_grid2, (3, 4, 5, 0, 1, 2))

    np.testing.assert_allclose(wave_grid1, wave_grid2, atol=1e-8)

test_wave()



def test_nabla_n_grid():
    pw_param = jr.pw.param_init(key, num_bands, 1, freq_mask)
    coeff = jr.pw.coeff(pw_param, freq_mask)
    occupation = jr.occupation.gamma(
      1, crystal.num_electron, num_bands=num_bands
    )

    r = jnp.array([0.1, 0.05, -0.6])
    nabla_density_r1 = jr.pw.nabla_density_grid(
      r, coeff, crystal.A, g_vecs, occupation
    )

    def nabla_n(r):
        def f(r):
            return jr.pw.density_r(
                r, coeff, crystal.A, g_vecs, occupation
            )
        return jax.grad(f)(r)

    nabla_density_r2 = nabla_n(r)

    np.testing.assert_allclose(nabla_density_r1, nabla_density_r2, atol=1e-8)

