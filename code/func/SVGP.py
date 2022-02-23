#!/usr/bin/env python

import numpy as np
import gpflow
import tensorflow as tf
from gpflow.ci_utils import ci_niter


class Matern25_aniso(gpflow.kernels.AnisotropicStationary):
    def K_d(self, d):
        sqrt5 = np.sqrt(5.0)
        d = tf.square(d)
        d = tf.reduce_sum(d, -1)
        d = tf.sqrt(d)
        tf.debugging.check_numerics(d, message='Checking distance matrix\n')

        return self.variance * (1.0 + sqrt5 * d + 5.0 / 3.0 *
                                tf.square(d)) * tf.exp(-sqrt5 * d)

def run_adam(model, iterations, tensor_data):
    """
    Utility function running the Adam optimizer

    :param model: GPflow model
    :param interations: number of iterations
    """
    # Create an Adam Optimizer action
    logf = []
    training_loss = model.training_loss_closure(tensor_data, compile=True)
    optimizer = tf.optimizers.Adam()

    @tf.function
    def optimization_step():
        optimizer.minimize(training_loss, model.trainable_variables)

    for step in range(iterations):
        print(model.trainable_variables)
        sys.stdout.flush()
        optimization_step()
        if step % 10 == 0:
            elbo = -training_loss().numpy()
            logf.append(elbo)
    return logf

def SVGP_wrap(locs, centers, thetaInit):
    tf.random.set_seed(123)
    n, d = locs.shape
    y = np.zeros((n, 1))
    theta = np.array(thetaInit)
    r = np.sqrt(theta[1: d + 1])
    l = 1 / r
    data = (locs, y)
    tensor_data = tuple(map(tf.convert_to_tensor, data))
    kernel = Matern25_aniso(variance=theta[0],
                            lengthscales=l)
    mdl = gpflow.models.SVGP(kernel, gpflow.likelihoods.Gaussian(),
                             centers, num_data=n, mean_function=None)
    gpflow.set_trainable(mdl.kernel, False)
    gpflow.set_trainable(mdl.inducing_variable, True)
    maxiter = ci_niter(1000)
    run_adam(mdl, maxiter, tensor_data)
    return np.array(mdl.inducing_variable.Z.numpy())
