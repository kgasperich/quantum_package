=====================
read_integral_complex
=====================

to use code with kpoints:

obtain and build as usual; install `read_integral_complex` and `Full_CI_ZMQ` plugins

``git clone https://github.com/kgasperich/quantum_package.git``

``cd quantum_package``

``./configure config/ifort.cfg``

``source quantum_package.rc``
``qp_module.py install read_integral_complex Full_CI_ZMQ``

``ninja``
build qp

install read_integral_complex and Full_CI_ZMQ modules

do pyscf calculation; call pyscf2QP to generate QP inputs

create ezfio with Gen_Ezfio_from_integral_complex.sh (this will also construct mo bielec map and save it to disk)

do cipsi with qp_run fci_zmq

Needed Modules
==============
.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.
Documentation
=============
.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.
