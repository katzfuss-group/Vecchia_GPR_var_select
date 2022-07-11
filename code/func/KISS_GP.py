import math
import torch
import gpytorch
import numpy


class GPRegressionModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(GPRegressionModel, self).__init__(train_x, train_y, likelihood)

        # SKI requires a grid size hyperparameter. This util can help with that
        # We're setting Kronecker structure to False because we're using an additive structure decomposition
        grid_size = int(gpytorch.utils.grid.choose_grid_size(train_x, kronecker_structure=False))

        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.AdditiveStructureKernel(
            gpytorch.kernels.ScaleKernel(
                gpytorch.kernels.GridInterpolationKernel(
                    gpytorch.kernels.MaternKernel(nu=2.5),
                    grid_size=grid_size, num_dims=1
                )
            ), num_dims=train_x.shape[-1]
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def KISS_GP_wrap(XTrn, yTrn, XTst, yTst):
    XTrnTorch = torch.from_numpy(XTrn).to(torch.float32)
    yTrnTorch = torch.from_numpy(numpy.array(yTrn)).to(torch.float32)
    XTstTorch = torch.from_numpy(XTst).to(torch.float32)
    yTstTorch = torch.from_numpy(numpy.array(yTst)).to(torch.float32)
    # Train
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = GPRegressionModel(XTrnTorch, yTrnTorch, likelihood)
    model.train()
    likelihood.train()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.1)
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
    training_iterations = 100
    for i in range(training_iterations):
        optimizer.zero_grad()
        output = model(XTrnTorch)
        loss = -mll(output, yTrnTorch)
        loss.backward()
        optimizer.step()
    # Predict
    model.eval()
    likelihood.eval()
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        yPred = likelihood(model(XTstTorch)).mean.view([XTstTorch.shape[0]])
        rmseScr = torch.sqrt(torch.square(yTstTorch - yPred).mean() / torch.var(yTstTorch))
    return dict({"yTstPred": yPred.tolist(), "rmseScr": rmseScr.item()})

## A test here
# n = 1000
# train_x = torch.rand([n, 4])
# train_y = torch.sin(train_x[:, 0] * (4 * math.pi) + torch.randn(n) * 0.2) + \
#     torch.sin(train_x[:, 1] * (2 * math.pi) + torch.randn(n) * 0.2) + \
#     torch.sin(train_x[:, 2] * (3 * math.pi) + torch.randn(n) * 0.2)
# nTst = 1000
# xTst = torch.rand([nTst, 4])
# yTst = torch.sin(xTst[:, 0] * (4 * math.pi) + torch.randn(n) * 0.2) + \
#     torch.sin(xTst[:, 1] * (2 * math.pi) + torch.randn(n) * 0.2) + \
#     torch.sin(xTst[:, 2] * (3 * math.pi) + torch.randn(n) * 0.2)
# KISS_GP_wrap(train_x.numpy(), train_y.numpy(), xTst.numpy(), yTst.numpy())
