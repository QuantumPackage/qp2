Selected Configuration Interaction
==================================

.. default-role:: cite

These methods rely on the same principle as the usual |CI| approaches, except
that determinants aren't chosen *a priori* based on an occupation or
excitation criterion, but selected *on the fly* among the entire set of
determinants based on their estimated contribution to the |FCI| wave function.
It has been noticed long ago that, even inside a predefined subspace of
determinants, only a small number significantly contributes to the wave
function. `Bytautas_2009,Anderson_2018` Therefore, an *on-the-fly*
selection of determinants is a rather natural idea that has been proposed
in the late 60's by Bender and Davidson `Bender_1969` as well as Whitten
and Hackmeyer. `Whitten_1969`

The approach we are using in the |qp| is based on |CIPSI| developed by Huron,
Rancurel and Malrieu, `Huron_1973` that iteratively selects *external*
determinants (determinants which are not present in the variational space)
using a perturbative criterion.

There is however a computational downside. In *a priori* selected
methods, the rule by which determinants are selected is known *a
priori*, and therefore, one can map a particular determinant to some row or
column index. `Knowles_1984` As a consequence, it can be systematically
determined to which matrix element of :math:`\hat H` a two-electron integral
contributes. This allows for the implementation of so-called
*integral-driven* methods, that work essentially by iterating over
integrals and are very fast.

On the contrary, in selected methods an explicit list of determinants has to be
kept, and there is no immediate way to know whether a determinant has been
selected, or what its index is in the list. Consequently, a
*determinant-driven* approach will be used, in which the loops run over
determinants rather than integrals. This can be a lot more computationally
expensive since the number of determinants is typically much larger than the
number of integrals.

What makes *determinant-driven* approaches possible here is:

- the fact that selected |CI| methods will keep the number of determinants small
  enough, orders of magnitude smaller than in *a priori* selected methods for
  wave functions with equal energies,
- an efficient way to compare determinants in order to extract the
  corresponding excitation operators `Scemama_2013`,
- an intense filtering of the internal space to avoid as much as possible
  determinant comparisons of disconnected determinants,
- a fast retrieval of the corresponding two-electron integrals in memory.


Simple Algorithm
----------------

.. default-role:: math 

.. |SetDI| replace:: `\{|D_I\rangle\}^{(n)}`
.. |Psi_n| replace:: `|\Psi^{(n)}\rangle`
.. |kalpha_star| replace:: `\{ |\alpha \rangle \}_\star ^{(n)}`
.. |ealpha| replace:: `e_\alpha`
.. |EPT| replace:: `E_\text{PT2}`

The variational wave function |Psi_n| is defined over a set of determinants
|SetDI| in which we diagonalize |H|.

.. math::

   |\Psi^{(n)}\rangle = \sum_{I} c_I^{(n)} |D_I\rangle


The determinants in |SetDI| will be characterized as **internal**.

#. For all **external** determinants |kalpha| `\notin` |SetDI|, compute the
   Epstein-Nesbet second-order perturbative contribution to the energy

   .. math::

      e_\alpha = \frac{ \langle \Psi^{(n)}| {\hat H} | \alpha \rangle^2 }{E^{(n)} - \langle \alpha | {\hat H} | \alpha \rangle }.

   `E^{(n)}` is the variational energy of the wave function at the current
   iteration. Note that another perturbation theory could be used to estimate
   |ealpha|.

#. An estimate of the total missing correlation energy can be computed
   by summing all the |ealpha| contributions

   .. math::

      E_\text{PT2} & = \sum_{\alpha} e_\alpha \\
      E_\text{FCI} & \approx E + E_\text{PT2} 

#. |kalpha_star|, the subset of determinants |kalpha| with the largest
   contributions |ealpha|, is added to the variational space

   .. math::

      \{ |D_I \rangle \}^{(n+1)} = \{|D_I\rangle\}^{(n)} \cup \{ |\alpha\rangle \}_\star^{(n)}


#. Go to iteration n+1, or exit on some criterion (number of determinants in
   the wave function, low |EPT|, ...).


The same procedure applies to any state, so both ground and excited states
can be treated.


Stochastic approximations for the selection and the computation of |EPT|
------------------------------------------------------------------------

The simple algorithm would be too slow to make calculations possible. Instead,
the |QP| uses the stochastic (hybrid) algorithm of :cite:`Garniron_2017.2` to
compute the |EPT| efficiently and to select on-the-fly the best Slater
determinants.

The selection step then introduces no extra cost with respect to the |EPT|
calculation, and the |EPT| itself is unbiased but carries a statistical error
bar that converges rapidly. In the hybrid algorithm, the stochastic part is
gradually reduced by increasing the deterministic part until a target
relative error is reached (:option:`perturbation pt2_relative_error`).


Deterministic approximations for the selection
----------------------------------------------

The following description corresponds to the fully deterministic algorithm,
used when :option:`perturbation do_pt2` is set to |false|. For small or
moderate wave functions it is cheaper than the hybrid algorithm, since every
generated candidate is evaluated exactly once, without the overhead of the
sampling machinery, and the hybrid algorithm itself resorts to this
deterministic evaluation only for the very smallest wave functions. The
|EPT| estimate is, however, biased by the truncation of the wave function,
so the extrapolated |FCI| energy is less accurate than with the hybrid
algorithm. It is suited to runs whose only goal is to generate a wave
function, for instance as a trial wave function for |QMC| calculations. It
also uses the notions of **generator** and **selector** determinants, which
are more general than the |CIPSI| algorithm (which targets the |FCI|) and
can be used to realize many kinds of selected |CI|.


We define **generator** determinants, as determinants of the internal space
from which the |kalpha| are generated.
We then define **selector** determinants, the sum of which forms a truncated
wave function used in the computation of |ealpha|.

For calculations in the |FCI| space, the determinants are sorted by decreasing
`|c_I|^2`, and thresholds are used on the squared norm of the wave function.
The default of :option:`determinants threshold_generators` is 0.999 for the
generators. The threshold for the selectors is not an independent parameter:
it is set to the square root of the generator threshold, which gives about
0.9995 by default.

This is the 3-class |CIPSI| approximation, :cite:`Evangelisti_1983`: instead
of generating all possible |kalpha|, only a subset likely to be selected is
generated.


The computation of |EPT| using a truncated wave function is biased,
so if an accurate estimate of the |FCI| energy is desired, it is preferable
to recompute |EPT| with the hybrid deterministic/stochastic algorithm
:cite:`Garniron_2017b` which is unbiased (this is the default).


Modifying the selection space
-----------------------------

By changing the definition of generators, and the rules for the generation of
the |kalpha|, it is easy to define selected variants of traditional |CI| methods.

For example, if one defines the |HF| determinant as the only generator,
one will produce a selected |CISD|. If one also changes the rules for the generation
to generate only the double excitations, one will have a selected |CID|.

The generators can also be chosen as determinants belonging to a |CAS|. If the
rules allow only for excitations inside the |CAS|, we obtain a selected
|CAS| |CI|. If the rules allow single and double excitations from the |CAS|
in the |FCI| space, we obtain a selected |CAS-SD|. And if one adds the rule
to prevent double excitations with two holes and two particles outside of the
active space, one obtains a selected |DDCI| method.

All such things can be done very easily when programming the |qp|.

.. bibliography::
   :style: unsrt
   :labelprefix: A


