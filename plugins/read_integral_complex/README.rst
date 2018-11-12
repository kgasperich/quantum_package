=====================
read_integral_complex
=====================

to use code with kpoints:

obtain and build as usual

install ``read_integral_complex`` and ``Full_CI_ZMQ`` plugins

``git clone https://github.com/kgasperich/quantum_package.git``

``cd quantum_package``

``./configure config/ifort.cfg`` (or other config file)

``source quantum_package.rc``

``qp_module.py install read_integral_complex Full_CI_ZMQ``

``ninja``


Perform pyscf calculation and use ``pyscf2QP`` function from ``MolPyscfToQPkpts.py`` to generate integrals and parameters. (see ``example.py`` in this directory)

generate ezfio:

``./Gen_Ezfio_from_integral_complex.sh name.ezfio``


do CIPSI:

``qp_run fci_zmq name.ezfio``


Needed Modules
==============
.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.
Documentation
=============
.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.
