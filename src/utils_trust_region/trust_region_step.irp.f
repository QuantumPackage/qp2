! Trust region

! *Compute the next step with the trust region algorithm*

! The Newton method is an iterative method to find a minimum of a given
! function. It uses a Taylor series truncated at the second order of the
! targeted function and gives its minimizer. The minimizer is taken as
! the new position and the same thing is done. And by doing so
! iteratively the method find a minimum, a local or global one depending
! of the starting point and the convexity/nonconvexity of the targeted
! function.  

! The goal of the trust region is to constrain the step size of the
! Newton method in a certain area around the actual position, where the 
! Taylor series is a good approximation of the targeted function. This
! area is called the "trust region".

! In addition, in function of the agreement between the Taylor
! development of the energy and the real energy, the size of the trust
! region will be updated at each iteration. By doing so, the step sizes
! are not too larges. In addition, since we add a criterion to cancel the
! step if the energy increases (more precisely if rho < 0.1), so it's
! impossible to diverge. \newline

! References: \newline
! Nocedal & Wright, Numerical Optimization, chapter 4 (1999), \newline
! https://link.springer.com/book/10.1007/978-0-387-40065-5, \newline
! ISBN: 978-0-387-40065-5 \newline

! By using the first and the second derivatives, the Newton method gives
! a step:
! \begin{align*}
!   \textbf{x}_{(k+1)}^{\text{Newton}} = - \textbf{H}_{(k)}^{-1} \cdot 
!   \textbf{g}_{(k)}
! \end{align*}
! which leads to the minimizer of the Taylor series.
! !!! Warning: the Newton method gives the minimizer if and only if
! $\textbf{H}$ is positive definite, else it leads to a saddle point !!!
! But we want a step $\textbf{x}_{(k+1)}$ with a constraint on its (euclidian) norm:
! \begin{align*}
!   ||\textbf{x}_{(k+1)}|| \leq \Delta_{(k+1)}
! \end{align*}
! which is equivalent to 
! \begin{align*}
!   \textbf{x}_{(k+1)}^T \cdot \textbf{x}_{(k+1)} \leq \Delta_{(k+1)}^2
! \end{align*}

! with: \newline
! $\textbf{x}_{(k+1)}$ is the step for the k+1-th iteration (vector of
! size n) \newline
! $\textbf{H}_{(k)}$ is the hessian at the k-th iteration (n by n
! matrix) \newline
! $\textbf{g}_{(k)}$ is the gradient at the k-th iteration (vector of
! size n) \newline
! $\Delta_{(k+1)}$ is the trust radius for the (k+1)-th iteration
! \newline

! Thus we want to constrain the step size $\textbf{x}_{(k+1)}$ into a
! hypersphere of radius $\Delta_{(k+1)}$.\newline

! So, if $||\textbf{x}_{(k+1)}^{\text{Newton}}|| \leq \Delta_{(k)}$ and
! $\textbf{H}$ is positive definite, the
! solution is the step given by the Newton method
! $\textbf{x}_{(k+1)} = \textbf{x}_{(k+1)}^{\text{Newton}}$.
! Else we have to constrain the step size. For simplicity we will remove
! the index $_{(k)}$ and $_{(k+1)}$. To restict the step size, we have
! to put a constraint on $\textbf{x}$ with a Lagrange multiplier.
! Starting from the Taylor series of a function E (here, the energy)
! truncated at the 2nd order, we have:
! \begin{align*}
!   E(\textbf{x}) =  E +\textbf{g}^T \cdot \textbf{x} + \frac{1}{2}
!   \cdot \textbf{x}^T \cdot \textbf{H} \cdot \textbf{x} +
!   \mathcal{O}(\textbf{x}^2)
! \end{align*} 

! With the constraint on the norm of $\textbf{x}$ we can write the
! Lagrangian
! \begin{align*}
!   \mathcal{L}(\textbf{x},\lambda) = E + \textbf{g}^T \cdot \textbf{x}
!   + \frac{1}{2} \cdot \textbf{x}^T \cdot \textbf{H} \cdot \textbf{x} 
!   + \frac{1}{2} \lambda (\textbf{x}^T \cdot \textbf{x} - \Delta^2)
! \end{align*}
! Where: \newline
! $\lambda$ is the Lagrange multiplier \newline
! $E$ is the energy at the k-th iteration $\Leftrightarrow
! E(\textbf{x} = \textbf{0})$ \newline

! To solve this equation, we search a stationary point where the first
! derivative of $\mathcal{L}$ with respect to $\textbf{x}$ becomes 0, i.e.
! \begin{align*}
!   \frac{\partial \mathcal{L}(\textbf{x},\lambda)}{\partial \textbf{x}}=0
! \end{align*}

! The derivative is:
! \begin{align*}
!   \frac{\partial \mathcal{L}(\textbf{x},\lambda)}{\partial \textbf{x}}
!   = \textbf{g} + \textbf{H} \cdot \textbf{x} + \lambda \cdot \textbf{x} 
! \end{align*}

! So, we search $\textbf{x}$ such as:
! \begin{align*}
! \frac{\partial \mathcal{L}(\textbf{x},\lambda)}{\partial \textbf{x}}
! = \textbf{g} + \textbf{H} \cdot \textbf{x} + \lambda \cdot \textbf{x} = 0
! \end{align*}

! We can rewrite that as:
! \begin{align*}
!   \textbf{g} + \textbf{H} \cdot \textbf{x} + \lambda \cdot \textbf{x} 
!   = \textbf{g} + (\textbf{H} +\textbf{I} \lambda) \cdot \textbf{x} = 0
! \end{align*}
! with $\textbf{I}$ is the identity matrix. 

! By doing so, the solution is:
! \begin{align*}
!   (\textbf{H} +\textbf{I} \lambda) \cdot \textbf{x}= -\textbf{g}
! \end{align*}
! \begin{align*}
!   \textbf{x}= - (\textbf{H} + \textbf{I} \lambda)^{-1} \cdot \textbf{g}
! \end{align*}
! with $\textbf{x}^T \textbf{x} = \Delta^2$.

! We have to solve this previous equation to find this $\textbf{x}$ in the
! trust region, i.e. $||\textbf{x}|| = \Delta$. Now, this problem is
! just a one dimension problem because we can express $\textbf{x}$ as a
! function of $\lambda$: 
! \begin{align*}
!   \textbf{x}(\lambda) = - (\textbf{H} + \textbf{I} \lambda)^{-1} \cdot \textbf{g}
! \end{align*}

! We start from the fact that the hessian is diagonalizable. So we have: 
! \begin{align*}
!   \textbf{H} = \textbf{W} \cdot \textbf{h} \cdot \textbf{W}^T
! \end{align*}
! with: \newline
! $\textbf{H}$, the hessian matrix \newline
! $\textbf{W}$, the matrix containing the eigenvectors \newline
! $\textbf{w}_i$, the i-th eigenvector, i.e. i-th column of $\textbf{W}$ \newline
! $\textbf{h}$, the matrix containing the eigenvalues in ascending order \newline
! $h_i$, the i-th eigenvalue in ascending order \newline

! Now we use the fact that adding a constant on the diagonal just shifts
! the eigenvalues:
! \begin{align*}
!   \textbf{H} + \textbf{I} \lambda = \textbf{W} \cdot (\textbf{h} 
!   +\textbf{I} \lambda) \cdot \textbf{W}^T
! \end{align*}

! By doing so we can express $\textbf{x}$ as a function of $\lambda$
! \begin{align*}
!   \textbf{x}(\lambda) = - \sum_{i=1}^n \frac{\textbf{w}_i^T \cdot 
!   \textbf{g}}{h_i + \lambda} \cdot \textbf{w}_i
! \end{align*}
! with $\lambda \neq - h_i$.

! An interesting thing in our case is the norm of $\textbf{x}$,
! because we want $||\textbf{x}|| = \Delta$. Due to the orthogonality of
! the eigenvectors $\left\{\textbf{w} \right\} _{i=1}^n$ we have:
! \begin{align*}
!   ||\textbf{x}(\lambda)||^2 = \sum_{i=1}^n \frac{(\textbf{w}_i^T \cdot 
!   \textbf{g})^2}{(h_i + \lambda)^2}
! \end{align*}

! So the $||\textbf{x}(\lambda)||^2$ is just a function of $\lambda$. 
! And if we study the properties of this function we see that: 
! \begin{align*}
!   \lim_{\lambda\to\infty} ||\textbf{x}(\lambda)|| = 0
! \end{align*}
! and if $\textbf{w}_i^T \cdot \textbf{g} \neq 0$: 
! \begin{align*}
!   \lim_{\lambda\to -h_i} ||\textbf{x}(\lambda)|| = + \infty
! \end{align*}

! From these limits and knowing that $h_1$ is the lowest eigenvalue, we
! can conclude that $||\textbf{x}(\lambda)||$ is a continuous and
! strictly decreasing function on the interval $\lambda \in
! (-h_1;\infty)$. Thus, there is one $\lambda$ in this interval which
! gives $||\textbf{x}(\lambda)|| = \Delta$, consequently there is one
! solution. 

! Since $\textbf{x} = - (\textbf{H} + \lambda \textbf{I})^{-1} \cdot
! \textbf{g}$ and we want to reduce the norm of $\textbf{x}$, clearly,
! $\lambda > 0$ ($\lambda = 0$ is the unconstraint solution). But the
! Newton method is only defined for a positive definite hessian matrix,
! so $(\textbf{H} + \textbf{I} \lambda)$ must be positive
! definite. Consequently, in the case where $\textbf{H}$ is not positive
! definite, to ensure the positive definiteness, $\lambda$ must be
! greater than $- h_1$.
! \begin{align*}
!   \lambda > 0 \quad \text{and} \quad \lambda \geq - h_1
! \end{align*} 

! From that there are five cases:
! - if $\textbf{H}$ is positive definite, $-h_1 < 0$, $\lambda \in (0,\infty)$
! - if $\textbf{H}$ is not positive definite and $\textbf{w}_1^T \cdot
!   \textbf{g} \neq 0$, $(\textbf{H} + \textbf{I}
!   \lambda)$ 
!   must be positve definite, $-h_1 > 0$, $\lambda \in (-h_1, \infty)$
! - if $\textbf{H}$ is not positive definite , $\textbf{w}_1^T \cdot
!   \textbf{g} = 0$ and $||\textbf{x}(-h_1)|| > \Delta$ by removing
!   $j=1$ in the sum, $(\textbf{H} + \textbf{I} \lambda)$ must be
!   positive definite, $-h_1 > 0$, $\lambda \in (-h_1, \infty$)
! - if $\textbf{H}$ is not positive definite , $\textbf{w}_1^T \cdot
!   \textbf{g} = 0$ and $||\textbf{x}(-h_1)|| \leq \Delta$ by removing
!   $j=1$ in the sum, $(\textbf{H} + \textbf{I} \lambda)$ must be
!   positive definite, $-h_1 > 0$, $\lambda = -h_1$). This case is
!   similar to the case where $\textbf{H}$ and $||\textbf{x}(\lambda =
!   0)|| \leq \Delta$
!   but we can also add to $\textbf{x}$, the first eigenvector $\textbf{W}_1$
!   time a constant to ensure the condition $||\textbf{x}(\lambda =
!   -h_1)|| = \Delta$ and escape from the saddle point

! Thus to find the solution, we can write: 
! \begin{align*}
!   ||\textbf{x}(\lambda)|| = \Delta
! \end{align*}
! \begin{align*}
!   ||\textbf{x}(\lambda)|| - \Delta = 0
! \end{align*}

! Taking the square of this equation
! \begin{align*}
!   (||\textbf{x}(\lambda)|| - \Delta)^2 = 0
! \end{align*}
! we have a function with one minimum for the optimal $\lambda$.
! Since we have the formula of $||\textbf{x}(\lambda)||^2$, we solve
! \begin{align*}
!   (||\textbf{x}(\lambda)||^2 - \Delta^2)^2 = 0
! \end{align*}

! But in practice, it is more effective to solve:
! \begin{align*}
!   (\frac{1}{||\textbf{x}(\lambda)||^2} - \frac{1}{\Delta^2})^2 = 0
! \end{align*}

! To do that, we just use the Newton method with "trust_newton" using
! first and second derivative of $(||\textbf{x}(\lambda)||^2 -
! \Delta^2)^2$ with respect to $\textbf{x}$.
! This will give the optimal $\lambda$ to compute the
! solution $\textbf{x}$ with the formula seen previously:
! \begin{align*}
!   \textbf{x}(\lambda) = - \sum_{i=1}^n \frac{\textbf{w}_i^T \cdot
!   \textbf{g}}{h_i + \lambda} \cdot \textbf{w}_i
! \end{align*}

! The solution $\textbf{x}(\lambda)$ with the optimal $\lambda$ is our
! step to go from the (k)-th to the (k+1)-th iteration, is noted $\textbf{x}^*$.




! Evolution of the trust region

! We initialize the trust region at the first iteration using a radius
! \begin{align*}
!   \Delta = ||\textbf{x}(\lambda=0)||
! \end{align*}

! And for the next iteration the trust region will evolves depending of
! the agreement of the energy prediction based on the Taylor series
! truncated at the 2nd order and the real energy. If the Taylor series
! truncated at the 2nd order represents correctly the energy landscape
! the trust region will be extent else it will be reduced. In order to
! mesure this agreement we use the ratio rho cf. "rho_model" and
! "trust_e_model". From that we use the following values:
! - if $\rho \geq 0.75$, then $\Delta = 2 \Delta$,
! - if $0.5 \geq \rho < 0.75$, then $\Delta = \Delta$, 
! - if $0.25 \geq \rho < 0.5$, then $\Delta = 0.5 \Delta$, 
! - if $\rho < 0.25$, then $\Delta = 0.25 \Delta$.

! In addition, if $\rho < 0.1$ the iteration is cancelled, so it
! restarts with a smaller trust region until the energy decreases.




! Summary

! To summarize, knowing the hessian (eigenvectors and eigenvalues), the
! gradient and the radius of the trust region we can compute the norm of
! the Newton step  
! \begin{align*}
!   ||\textbf{x}(\lambda = 0)||^2 = ||- \textbf{H}^{-1} \cdot \textbf{g}||^2 = \sum_{i=1}^n 
!   \frac{(\textbf{w}_i^T \cdot \textbf{g})^2}{(h_i + \lambda)^2}, \quad h_i \neq 0
! \end{align*}

! - if $h_1 \geq 0$, $||\textbf{x}(\lambda = 0)|| \leq \Delta$ and
!   $\textbf{x}(\lambda=0)$ is in the trust region and it is not
!   necessary to put a constraint on $\textbf{x}$, the solution is the
!   unconstrained one, $\textbf{x}^* = \textbf{x}(\lambda = 0)$.
! - else if $h_1 < 0$, $\textbf{w}_1^T \cdot \textbf{g} = 0$ and
!   $||\textbf{x}(\lambda = -h_1)|| \leq \Delta$ (by removing $j=1$ in
!   the sum), the solution is $\textbf{x}^* = \textbf{x}(\lambda =
!   -h_1)$, similarly to the previous case.
!   But we can add to $\textbf{x}$, the first eigenvector $\textbf{W}_1$
!   time a constant to ensure the condition $||\textbf{x}(\lambda =
!   -h_1)|| = \Delta$ and escape from the saddle point
! - else if $h_1 < 0$ and $\textbf{w}_1^T \cdot \textbf{g} \neq 0$ we
!   have to search $\lambda \in (-h_1, \infty)$ such as
!   $\textbf{x}(\lambda) = \Delta$ by solving with the Newton method 
!   \begin{align*}
!     (||\textbf{x}(\lambda)||^2 - \Delta^2)^2 = 0
!   \end{align*}
!   or
!   \begin{align*}
!     (\frac{1}{||\textbf{x}(\lambda)||^2} - \frac{1}{\Delta^2})^2 = 0
!   \end{align*}
!   which is numerically more stable. And finally compute 
!   \begin{align*}
!     \textbf{x}^* = \textbf{x}(\lambda) = - \sum_{i=1}^n \frac{\textbf{w}_i^T \cdot
!     \textbf{g}}{h_i + \lambda} \cdot \textbf{w}_i
!   \end{align*}
! - else if $h_1 \geq 0$ and $||\textbf{x}(\lambda = 0)|| > \Delta$ we
!   do exactly the same thing that the previous case but we search
!   $\lambda \in (0, \infty)$ 
! - else if $h_1 < 0$ and $\textbf{w}_1^T \cdot \textbf{g} = 0$ and
!   $||\textbf{x}(\lambda = -h_1)|| > \Delta$ (by removing $j=1$ in the
!   sum), again we do exactly the same thing that the previous case
!   searching $\lambda \in (-h_1, \infty)$.
  

! For the cases where $\textbf{w}_1^T \cdot \textbf{g} = 0$ it is not
! necessary in fact to remove the $j = 1$ in the sum since the term
! where $h_i - \lambda < 10^{-6}$ are not computed.

! After that, we take this vector $\textbf{x}^*$, called "x", and we do
! the transformation to an antisymmetric matrix $\textbf{X}$, called
! m_x. This matrix $\textbf{X}$ will be used to compute a rotation
! matrix $\textbf{R}= \exp(\textbf{X})$ in "rotation_matrix".

! NB: 
! An improvement can be done using a elleptical trust region.




! Code

! Provided:
! | mo_num | integer | number of MOs |

! Cf. qp_edit in orbital optimization section, for some constants/thresholds

! Input:
! | m         | integer          | number of MOs                                |
! | n         | integer          | m*(m-1)/2                                       |
! | H(n, n)   | double precision | hessian                                         |
! | v_grad(n) | double precision | gradient                                        |
! | e_val(n)  | double precision | eigenvalues of the hessian                      |
! | W(n, n)   | double precision | eigenvectors of the hessian                     |
! | rho       | double precision | agreement between the model and the reality,    |
! |           |                  | represents the quality of the energy prediction |
! | nb_iter   | integer          | number of iteration                             |

! Input/Ouput:
! | delta | double precision | radius of the trust region |

! Output:
! | x(n)      | double precision | vector containing the step |

! Internal:
! | accu          | double precision | temporary variable to compute the step       |
! | lambda        | double precision | lagrange multiplier                          |
! | trust_radius2 | double precision | square of the radius of the trust region     |
! | norm2_x       | double precision | norm^2 of the vector x                       |
! | norm2_g       | double precision | norm^2 of the vector containing the gradient |
! | tmp_wtg(n)    | double precision | tmp_wtg(i) = w_i^T . g                       |
! | i, j, k       | integer          | indexes                                      |

! Function:
! | dnrm2                   | double precision | Blas function computing the norm       |
! | f_norm_trust_region_omp | double precision | compute the value of norm(x(lambda)^2) |


subroutine trust_region_step(n,nb_iter,v_grad,rho,e_val,w,x,delta)

  include 'pi.h'

  BEGIN_DOC
  ! Compuet the step in the trust region
  END_DOC

  implicit none

  ! Variables

  ! in
  integer, intent(in)             :: n
  double precision, intent(in)    :: v_grad(n), rho
  integer, intent(inout)          :: nb_iter
  double precision, intent(in)    :: e_val(n), w(n,n)

  ! inout
  double precision, intent(inout) :: delta

  ! out
  double precision, intent(out)   :: x(n)

  ! Internal
  double precision                :: accu, lambda, trust_radius2
  double precision                :: norm2_x, norm2_g
  double precision, allocatable   :: tmp_wtg(:)
  integer                         :: i,j,k
  double precision                :: t1,t2,t3
  integer                         :: n_neg_eval


  ! Functions
  double precision                :: ddot, dnrm2
  double precision                :: f_norm_trust_region_omp

  print*,''
  print*,'=================='
  print*,'---Trust_region---'
  print*,'=================='

  call wall_time(t1)

  ! Allocation
  allocate(tmp_wtg(n))

! Initialization and norm

! The norm of the step size will be useful for the trust region
! algorithm. We start from a first guess and the radius of the trust
! region will evolve during the optimization.

! avoid_saddle is actually a test to avoid saddle points


! Initialization of the Lagrange multiplier
lambda = 0d0

! List of w^T.g, to avoid the recomputation
tmp_wtg = 0d0
do j = 1, n
  do i = 1, n
    tmp_wtg(j) = tmp_wtg(j) + w(i,j) * v_grad(i)
  enddo
enddo

! Replacement of the small tmp_wtg corresponding to a negative eigenvalue
! in the case of avoid_saddle
if (avoid_saddle .and. e_val(1) < - thresh_eig) then
  i = 2
  ! Number of negative eigenvalues
  do while (e_val(i) < - thresh_eig)
    if (tmp_wtg(i) < thresh_wtg2) then
      if (version_avoid_saddle == 1) then  
        tmp_wtg(i) = 1d0
      elseif  (version_avoid_saddle == 2) then  
        tmp_wtg(i) = DABS(e_val(i))
      elseif  (version_avoid_saddle == 3) then  
        tmp_wtg(i) = dsqrt(DABS(e_val(i)))
      else
        tmp_wtg(i) = thresh_wtg2
      endif
    endif
    i = i + 1
  enddo

  ! For the fist one it's a little bit different
  if (tmp_wtg(1) < thresh_wtg2) then 
    tmp_wtg(1) = 0d0
  endif

endif 

! Norm^2 of x, ||x||^2
norm2_x = f_norm_trust_region_omp(n,e_val,tmp_wtg,0d0)
! We just use this norm for the nb_iter = 0 in order to initialize the trust radius delta
! We don't care about the sign of the eigenvalue we just want the size of the step in a normal Newton-Raphson algorithm
! Anyway if the step is too big it will be reduced
print*,'||x||^2 :', norm2_x

! Norm^2 of the gradient, ||v_grad||^2
norm2_g = (dnrm2(n,v_grad,1))**2
print*,'||grad||^2 :', norm2_g

! Trust radius initialization

!     At the first iteration (nb_iter = 0) we initialize the trust region
!     with the norm of the step generate by the Newton's method ($\textbf{x}_1 =
!     (\textbf{H}_0)^{-1} \cdot \textbf{g}_0$,
!     we compute this norm using f_norm_trust_region_omp as explain just
!     below) 


! trust radius
if (nb_iter == 0) then
   trust_radius2 = norm2_x 
   ! To avoid infinite loop of cancellation of this first step
   ! without changing delta
   nb_iter = 1

   ! Compute delta, delta = sqrt(trust_radius)
   delta = dsqrt(trust_radius2)
endif

! Modification of the trust radius

! In function of rho (which represents the agreement between the model
! and the reality, cf. rho_model) the trust region evolves. We update
! delta (the radius of the trust region).

! To avoid too big trust region we put a maximum size.


! Modification of the trust radius in function of rho
if (rho >= 0.75d0) then
   delta = 2d0 * delta
elseif (rho >= 0.5d0) then
   delta = delta
elseif (rho >= 0.25d0) then
   delta = 0.5d0 * delta
else
   delta = 0.25d0 * delta
endif

! Maximum size of the trust region
!if (delta > 0.5d0 * n * pi) then
!  delta = 0.5d0 * n * pi
!  print*,'Delta > delta_max, delta = 0.5d0 * n * pi'
!endif

if (delta > 1d10) then
  delta = 1d10
endif

print*, 'Delta :', delta

! Calculation of the optimal lambda

! We search the solution of $(||x||^2 - \Delta^2)^2 = 0$
! - If $||\textbf{x}|| > \Delta$  or $h_1 < 0$ we have to add a constant
!   $\lambda > 0 \quad \text{and} \quad \lambda > -h_1$
! - If $||\textbf{x}|| \leq \Delta$ and $h_1 \geq 0$ the solution is the
!   unconstrained one, $\lambda = 0$

! You will find more details at the beginning


! By giving delta, we search (||x||^2 - delta^2)^2 = 0
! and not (||x||^2 - delta)^2 = 0

! Research of lambda to solve ||x(lambda)|| = Delta 

! Display
print*, 'e_val(1) = ', e_val(1)
print*, 'w_1^T.g =', tmp_wtg(1)

! H positive definite 
if (e_val(1) > - thresh_eig) then
  norm2_x = f_norm_trust_region_omp(n,e_val,tmp_wtg,0d0)
  print*, '||x(0)||=', dsqrt(norm2_x)
  print*, 'Delta=', delta

  ! H positive definite, ||x(lambda = 0)|| <= Delta
  if (dsqrt(norm2_x) <= delta) then 
    print*, 'H positive definite, ||x(lambda = 0)|| <= Delta'
    print*, 'lambda = 0, no lambda optimization'
    lambda = 0d0

  ! H positive definite, ||x(lambda = 0)|| > Delta
  else
    ! Constraint solution
    print*, 'H positive definite, ||x(lambda = 0)|| > Delta' 
    print*,'Computation of the optimal lambda...'
    call trust_region_optimal_lambda(n,e_val,tmp_wtg,delta,lambda)
  endif

! H indefinite
else
  if (DABS(tmp_wtg(1)) < thresh_wtg) then
    norm2_x = f_norm_trust_region_omp(n,e_val,tmp_wtg, - e_val(1))
    print*, 'w_1^T.g <', thresh_wtg,', ||x(lambda = -e_val(1))|| =', dsqrt(norm2_x) 
  endif

  ! H indefinite, w_1^T.g = 0, ||x(lambda = -e_val(1))|| <= Delta 
  if (dsqrt(norm2_x) <= delta .and. DABS(tmp_wtg(1)) < thresh_wtg) then
    ! Add e_val(1) in order to have (H - e_val(1) I) positive definite
    print*, 'H indefinite, w_1^T.g = 0, ||x(lambda = -e_val(1))|| <= Delta'
    print*, 'lambda = -e_val(1), no lambda optimization'
    lambda = - e_val(1)

  ! H indefinite, w_1^T.g = 0, ||x(lambda = -e_val(1))|| > Delta
  ! and
  ! H indefinite, w_1^T.g =/= 0
  else
    ! Constraint solution/ add lambda
    if (DABS(tmp_wtg(1)) < thresh_wtg) then
       print*, 'H indefinite, w_1^T.g = 0, ||x(lambda = -e_val(1))|| > Delta'
    else
       print*, 'H indefinite, w_1^T.g =/= 0'
    endif
    print*, 'Computation of the optimal lambda...'
    call trust_region_optimal_lambda(n,e_val,tmp_wtg,delta,lambda)
    endif

endif

! Recomputation of the norm^2 of the step x
norm2_x = f_norm_trust_region_omp(n,e_val,tmp_wtg,lambda)
print*,''
print*,'Summary after the trust region:'
print*,'lambda:', lambda
print*,'||x||:', dsqrt(norm2_x)
print*,'delta:', delta

! Calculation of the step x

! x refers to $\textbf{x}^*$
! We compute x in function of lambda using its formula :
! \begin{align*}
! \textbf{x}^* = \textbf{x}(\lambda) = - \sum_{i=1}^n \frac{\textbf{w}_i^T \cdot \textbf{g}}{h_i 
! + \lambda} \cdot \textbf{w}_i
! \end{align*}


! Initialisation
x = 0d0

! Calculation of the step x

! Normal version
if (.not. absolute_eig) then

  do i = 1, n 
    if (DABS(e_val(i)) > thresh_eig .and. DABS(e_val(i)+lambda) > thresh_eig) then
      do j = 1, n
        x(j) = x(j) - tmp_wtg(i) * W(j,i) / (e_val(i) + lambda)
      enddo
    endif
  enddo

! Version to use the absolute value of the eigenvalues
else

  do i = 1, n 
    if (DABS(e_val(i)) > thresh_eig) then
      do j = 1, n
        x(j) = x(j) - tmp_wtg(i) * W(j,i) / (DABS(e_val(i)) + lambda)
      enddo
    endif
  enddo

endif

double precision :: beta, norm_x

! Test
! If w_1^T.g = 0, the lim of ||x(lambda)|| when lambda tend to -e_val(1)
! is not + infinity. So ||x(lambda=-e_val(1))|| < delta, we add the first 
! eigenvectors multiply by a constant to ensure the condition
! ||x(lambda=-e_val(1))|| = delta and escape the saddle point
if (avoid_saddle .and. e_val(1) < - thresh_eig) then
  if (tmp_wtg(1) < 1d-15 .and. (1d0 - dsqrt(norm2_x)/delta) > 1d-3 ) then

    ! norm of x
    norm_x = dnrm2(n,x,1)

    ! Computes the coefficient for the w_1 
    beta = delta**2 - norm_x**2

    ! Updates the step x
    x = x + W(:,1) * dsqrt(beta)

    ! Recomputes the norm to check
    norm_x = dnrm2(n,x,1)

    print*, 'Add w_1 * dsqrt(delta^2 - ||x||^2):'
    print*, '||x||', norm_x 
  endif
endif

! Transformation of x

! x is a vector of size n, so it can be write as a m by m
! antisymmetric matrix m_x cf. "mat_to_vec_index" and "vec_to_mat_index".


!  ! Step transformation vector -> matrix
!  ! Vector with n element -> mo_num by mo_num matrix
!  do j = 1, m
!     do i = 1, m
!        if (i>j) then
!           call mat_to_vec_index(i,j,k)
!           m_x(i,j) = x(k)
!        else
!           m_x(i,j) = 0d0
!        endif
!     enddo
!  enddo
!
!  ! Antisymmetrization of the previous matrix
!  do j = 1, m
!     do i = 1, m
!        if (i<j) then
!           m_x(i,j) = - m_x(j,i)
!        endif
!     enddo
!  enddo

! Deallocation, end


deallocate(tmp_wtg)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in trust_region:', t3
  print*,'======================'
  print*,'---End trust_region---'
  print*,'======================'
  print*,''

end
