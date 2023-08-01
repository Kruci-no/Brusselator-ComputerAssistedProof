
# Brusselator-ComputerAssistedProof

The code is used for computer-assisted proof of the existence of a periodic orbit and its attraction for the non-autonomous Chafee-Infante system:

```math
\begin{cases}
 u_t = d_1  u_{xx}- (B+1)u +u^2v + A \sin(x)   \text{ for } (x,t)\in (0,\pi)\times(0,\infty),\\
 v_t = d_2 v_{xx} + Bu - u^2v    \text{ for } (x,t)\in (0,\pi)\times(0,\infty),\\
u(t,x) = v(t,x )= 0 \text{ for }  (x,t)\in {0,\pi} \times (0,\infty),\\
u(0,x) = u^0(x), v(0,x) = v^0(x)\ \textrm{for}\ x\in(0,\pi).
\end{cases}
```

The system models autocatalytic chemical reactions. The numbers $A,B,d_1,d_2>0$ are the parameters of the system.

For the Brusselator system, as the solutions $u(t,x)$ and $v(t,x)$ satisfy the Dirichlet boundary conditions on $[0,\pi]$, they are represented as sums of sine Fourier series:

```math
u(t,x) = \sum_{i=1}^\infty u_i(t)\sin(i x),\quad
v(t,x) = \sum_{i=1}^\infty v_i(t)\sin(i x).
```

The code contains 3 programs:
- `findPeriodicPoint`,
- `makeSubsection`,
- `testSolver`.

These programs use data from files:
- `initialValue.txt` - It contains coefficients for sine odd series used as initial data for the programs. For example,

```
{3,4,0,-1,0.5,0.5}   
```

represents initial data in the form:

```math
u^0 = 3\sin(x) - 0\sin(3x) - 0.5\sin(5x),
v^0 = 4\sin(x) - 1\sin(3x) + 0.5\sin(5x).
```

The same rule of indexing applies to all vectors corresponding to the finite part of points from the phase space.
- Parameters ${d_1,d_2,B,A}$ are taken from file `params.txt` in the form:

```
{0.2, 0.02, 2, 1}
-------------
{d_1,d_2,B,A}
```

## sampleDyn.cpp

Program for numerically integrating the Brusselator system. 

- Integrating in the space of odd coefficients.
- The initial data is taken from the content of file `initialValue.txt`. The initial data determines the size of Gallerkin projection.
- By using the program:

```
sampleDyn 100 10 1000
```

where:
- `100` is the starting time,
- `10` is the duration of simulation, and
- `1000` is the number of steps,

the program will output $1000$ uniformly sampled points at times between $[100,110]$ starting from initial data and save them to file `output2.txt`.

## makeSubsection.cpp

Program for finding periodic orbits of the Brusselator system. It tries to find a fixed periodic point by iterating the Galerkin approximation of some Poincare map. Additionally, after finding the point, it proposes the section (along with intermediate sections if set).

- The initial data for searching is taken from the content of file `initialValue.txt`. This initial data determines the size of the Galerkin projection.
- In file `section.txt`, we set the options for this program:
```
0 0.7
1
1e-5
------
numVariable Value
numberOfSection
MaxAbsEigenValue
```
Here,
- $numVariable$ and $Value$ define the starting Poincare map and section for which we are searching for a fixed point. The $numVariable$ defines the number of variables that are set to be equal to $Value$ at this section.
- The $numberOfSection$ defines how many sections will be outputted by the program.
- The $MaxAbsEigenValue$ defines a threshold for eigenvalues, based on which corresponding eigenvectors are taken as coordinates on an outputted section. For other points, a different method is used.

The program also outputs subsections in 'subsections.txt', where origins are determined by trajectories obtained from numerical approximations of periodic orbits.

 ## CAProof.cpp

This program is designed to prove the existence of a periodic orbit. It does so by checking whether a defined set $X^0$ satisfies the following condition:

$$ P(X^0) \subset X^0,$$

where $P$ represents a Poincare map. If this condition is satisfied, it confirms the existence of a periodic orbit. The program utilizes a rigorous C0 algorithm for integrating partial differential equations (PDEs) to compute the image.

Here are some details about the program and its inputs:

- The section which defines the Poincare map is read from `subsections.txt`.
- In file `initialValuePoincare.txt`, we set:
```
{[-1e-4,1e-4],[-1e-7,1e-7]} {[-10,10],[4,4]}
-------------
L  {[C^-,C^+], [s,s]}

```
These values are used in defining $X^0$.
  
The set $X_0$ is defined as follows:
   $$X^0 = X_P^0 + X_Q^0,$$

where
$$X_P^0 = x_0+ Aq_0.$$
 - Parameter $L$ defines entries of $q_0$: first entry is 0 and the next ones are taken from $L$,

$$
    X_Q^0 = \left(\sum_{i=n+1}^\infty \frac{[C^- , C ^+ ]}{(2i-1)^s}\sin((2i-1)x),\sum_{i=n+1}^\infty \frac{[C^- , C ^+ ]}{(2i-1)^s}\sin((2i - 1)x)\right),
$$

Here, $n$ depends on the size of the sections in the `subsection.txt` file, and $C$ and $s$ are obtained from the file `initialValuePoincare.txt`.

- The $x_0$ comes from orgin of Poincare section.
- The number of modes explicitly represented.

The rigorous C0 integration is performed in the odd subspace of Fourier coefficients.

## Code Information

- The programs use the [CAPD library](http://capd.ii.uj.edu.pl/index.php) - a tool for nonrigorous and validated numerics for dynamical systems.
  - Used version of the library: 5.1.2.
  - The Makefile assumes that the CAPD library is located in the following position relative to the main directory:
   ```
   # directory where capd scripts are (e.g. capd-config)
   CAPDBINDIR = ../capd-capdDynSys-5.1.2/bin/
   
   ```
  
- File Algebra contains the structure infinite series used in the algorithm and implementation of operations on them.
- File VectorFieldMaker contains additional methods that allow producing vector field fields for Gallerkin projection in string form, which is used in CAPD IMap class.
- File SolverPDE.h contains methods and structs for rigorous C0 integration.
- File SetInOut.h mainly contains Input/Output functions.
- File DSubsection.h is used to save subsections.
- File BrusellatorTools.h contains functions used to set vector fields which are utilized in numerical and rigorous integration of Burusselator vectors.



