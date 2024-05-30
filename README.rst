.. image:: doc/logo.svg
  :alt: mc2err

=======================

``mc2err`` is an incomplete open-source C library for computing standardized error bars for Markov-chain data.
After three iterations of design and implementation, further development is on hiatus because I have other research priorities right now.

So far, this has been an interesting exercise in trying to balance the flexibility of data ingestion with the flexibility of statistical analysis.
In its current form, the library can ingest data from multiple realizations of a Markov chain asynchronously, and it allows for missing or conditional data.
I never settled on an approach to statistical analysis, and kept oscillating between hypothesis testing and maximum likelihood estimation.
With hypothesis testing, the outputs would be raw averages and covariances with equilibration times and autocovariance cutoffs determined by tests and user-specified p values.
With maximum likelihood estimation, an asymptotically stationary, banded covariance matrix would be fit to the data, and averages and covariances would be inferred from the 
model.
In either case, no problem-specific model would be assumed, and both forms of analysis would rely on asymptotic normality assumptions.

While I had mostly settled on hypothesis testing in the previous design iteration, I became increasingly dissatisfied by its interplay with the flexible data input.
It would be very possible for users to input data in a very inefficient way, such as one very long Markov chain and then a very large number of very short chains.
In such an example, the analysis of the long Markov chain might cause the equilibration time to be set to a point beyond the end of the short chains, causing none of their data 
to be used in the outputs.
In contrast, the short chains would still contribute to maximum likelihood estimation by refining the model for the short-time transient behavior of the Markov chain.
Ultimately, I think this flexible input design naturally pairs with maximum likelihood estimation, while hypothesis testing naturally pairs with a more rigid input scheme.
In the more rigid case, the analysis would simply ingest instantaneous averages from a pool of samples at same time for different realizations of the same Markov chain.
Such an accumulator would be much simpler than this library, and it is something I will consider separately from this project, although I might off both approaches in the same 
library eventually.
