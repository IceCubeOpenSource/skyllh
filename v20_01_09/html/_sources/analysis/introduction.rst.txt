.. _analysis_introduction:

************
Introduction
************

The user can find pre-defined IceCube log-likelihood analyses in "i3skyllh"
project.


Setup
=====

To set-up and run an analysis the following procedure applies:

1. Create an analysis instance.
2. Add the datasets and their PDF ratio instances via the
   :meth:`~.Analysis.add_dataset` method.
3. Construct the log-likelihood ratio function via the
   :meth:`~.Analysis.construct_llhratio` method.
4. Call the :meth:`.do_trial` or :meth:`.unblind` method to perform a
   random trial or to unblind the data. Both methods will fit the global
   fit parameters using the set up data. Finally, the test statistic
   is calculated via the :meth:`.calculate_test_statistic` method.
