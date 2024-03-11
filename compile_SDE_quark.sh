#!/bin/bash
cd ..
cmake SDE_quark
make

mv main_sde_pmts ./SDE_quark/
mv main_sde_iter ./SDE_quark/
mv main_complex_prop ./SDE_quark/
mv sde_iter_solver.cpython-* ./SDE_quark/python/
mv sde_ker.cpython-* ./SDE_quark/python/
mv ker_itg.cpython-* ./SDE_quark/python/

cd SDE_quark
