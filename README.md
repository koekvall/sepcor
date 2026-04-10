# sepcor

R package for likelihood-based inference with separable correlation matrices. The model is

Sigma = D (C2 x C1) D,

where C1 and C2 are correlation matrices, D is a diagonal matrix with positive entries, and x denotes the Kronecker product. The main functions are:

- `sepcor`: fits the model by maximum likelihood via block-coordinate descent
- `sepcor_se`: computes standard errors from the expected Fisher information

## Authors

Karl Oskar Ekvall (k.o.ekvall at gmail dot com)

## License

This project is licensed under the MIT License.
