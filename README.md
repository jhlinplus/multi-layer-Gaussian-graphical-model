# Penalized MLE of Multi-Layered Gaussian Graphical Models

This is the official code repository for paper titled **Penalized Maximum Likelihood Estimation of Multi-Layered Gaussian Graphical Models**, published in the _Journal of Machine Learning Research, 2016_. http://www.jmlr.org/papers/volume17/16-004/16-004.pdf

To cite this work: 
```
@article{Lin2016Penalized,
  author  = {Jiahe Lin and Sumanta Basu and Moulinath Banerjee and George Michailidis},
  title   = {Penalized Maximum Likelihood Estimation of Multi-layered Gaussian Graphical Models},
  journal = {Journal of Machine Learning Research},
  year    = {2016},
  volume  = {17},
  number  = {146},
  pages   = {1--51},
  url     = {http://jmlr.org/papers/v17/16-004.html}
}
```

For the time being, we provide `R` implementation of the proposed methodology, which is the one used during model dev. See `example.R` for a demo. 

### Notes

I personally found that R env has become a bit difficult to use from a maintainance standpoint (e.g., pkgs no longer being supported by a newer R version), and I plan to roll out a `Python` version of the algorithm when I find the time. 