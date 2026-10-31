.. _print_mos: 
 
.. program:: print_mos 
 
========= 
print_mos 
========= 
 
 
 
 
 Numerically integrates g*gtt*x**2 and g**2*x**2 over x = 0..5 
 (alpha = 2, 10000 steps); g = exp(-alpha*x**2) and gtt = g* 
 (4*alpha**2*x**2 - 4*alpha). Writes (x, g, gtt) to the unopened 
 unit 33, scales both sums by 4*pi, prints accu_norm and 
 -0.5*accu/accu_norm, plus ao_kinetic_integrals(1,1). 
 TODO(doc): no MOs are printed; what is the purpose of this program? 
 
 Needs: 
 
 .. hlist:: 
    :columns: 3 
 
    * :c:data:`ao_kinetic_integrals` 
