.. _unit_tests:

**********
Unit tests
**********

>>> import unittest

.. code-block:: python

    import unittest     

    class SelfDrivingCarTest(TestCase):
        def setUp(self):
            self.car = SelfDrivingCar()

To run all test we can use following command::

    python -m unittest discover