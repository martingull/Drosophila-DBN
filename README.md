# Drosophila-DBN
MATLAB codes of the Drosophila Dynamic Bayes Network (DBN) project

# Introduction
Being able to infer gene regulatory networks from spatio-temporal expression
data is a major problem in biology. This thesis proposes a new dynamic Bayes
networks approach, which we benchmark by using the well researched gap gene
problem of the Drosophila melanogaster, with the capability of realistically
inferring gene regulatory networks and producing high quality simulations. The
thesis solves practical issues, currently associated with spatio-temporal gene
inference, such as computational time and parameter fragility, while obtaining
a similar gene regulatory network and matrix as our ground truth network. The
proposed modelling framework computes the gene regulatory network in 10-15
second on a modern laptop. Effectively removing the computational barrier of
the problem and allowing for future gene regulatory networks of greater gene
count to be processed. Besides producing a gene regulatory matrix our method
also produces high quality simulations of the gene activation levels of the gap
gene problem. In addition, unlike many competing problem formulations, the
proposed model is probabilistic in nature, hence allowing statistical inference
to be made. Finally, using Bayesian statistics, we perform robustness tests on
the topology of our proposed gene regulatory network and our regulatory
weights.
