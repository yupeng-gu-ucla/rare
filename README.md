# The Co-Evolution Model for Social Network Evolving and Opinion Migration

Implementation of the paper [**RaRE: Social Rank Regulated Large-scale Network Embedding](http://yupenggu.me/papers/WWW18_RaRE.pdf)

##### Usage: 

```
java Main <options>
```

or 

```
java Main -help
```
to show the help command below.

##### options:

- -k: dimension of embedding; 2 by default
- -data: the location of data file
- -weighted: whether the network is weighted (0: unweighted; other: weighted); 0 by default
- -out: the directory of output files; "." (current directory) by default
- -lr: initial stepsize; 0.05 by default
- -lambda_r: coefficient of dr (first component), must be positive; 1.0 by default
- -lambda_z: coefficient of dz (second component), must be positive; 1.0 by default
- -lambda_0: bias term (third component), should be negative for sparse networks; -1.0 by default
- -alpha_r: parameter for power law prior (on r), recommended value: 1.0~2.0; 1.5 by default
- -reg_z: l2 regularization coefficient (on z); 1e-6 by default
- -tp: fraction of nodes for training, must be between 0 and 1; 0.9 by default
- -tolerance: stopping criterion on log-likelihood improvement (negative number for no tolerance); 1e-6 by default
- -nsw: negative sample weight, must be positive; 5 by default
- -iter: value of maximum edges to be sampled (in thousands); 10000 thousand by default
- -verbose: whether see verbose output or not (0: show limited outputs; other: show all outputs); 1 by default
- -h or -help: show this help command

The -data command is necessary.

### Citation:
```
@inproceedings{gu2018rare,
  title={RaRE: Social Rank Regulated Large-scale Network Embedding},
  author={Gu, Yupeng and Sun, Yizhou and Li, Yanen and Yang, Yang},
  booktitle={Proceedings of the 2018 World Wide Web Conference on World Wide Web},
  pages={359--368},
  year={2018},
  organization={International World Wide Web Conferences Steering Committee}
}
```

### Contact:
[Yupeng Gu](http://web.cs.ucla.edu/~ypgu/)

