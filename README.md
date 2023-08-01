# Brusselator-ComputerAssistedProof

The code is used for computer-assisted proof of the existence of a periodic orbit and its attraction for the non-autonomous Chafee-Infante system:

$$\begin{cases}
 u_t = d_1  u_{xx}- (B+1)u +u^2v + A \sin(x)   \text{ for } (x,t)\in (0,\pi)\times(0,\infty),\\
 v_t = d_2 v_{xx} + Bu - u^2v    \text{ for } (x,t)\in (0,\pi)\times(0,\infty),\\
u(t,x) = v(t,x )= 0 \text{ for }  (x,t)\in {0,\pi} \times (0,\infty),\\
u(0,x) = u^0(x), v(0,x) = v^0(x)\ \textrm{for}\ x\in(0,\pi).
\end{cases},$$

The system models autocatalytic chemical reactions. The  $A,B,d_1,d_2>0$ are the parameters of the system.

For the Brusselator system, as the solutions $u(t,x)$ and $v(t,x)$ satisfy the Dirichlet boundary conditions on $[0,\pi]$, they are represented as the sums of sine Fourier series

$$
u(t,x) = \sum_{i=1}^\infty u_i(t)\sin(i x),\quad
v(t,x) = \sum_{i=1}^\infty v_i(t)\sin(i x).
$$

The code contains 3 programs:

- `findPeriodicPoint`,
- `makeSubsection`,
- `testSolver`.

The programs use data from files:
- `initialValue.txt` - It contains coefficients for sine odd series used as initial data for the programs. For example,

```
{3,4,0,-1,0.5,0.5}   
```

represents initial data in the form

$$
u^0 = 3\sin(x) - 0sin(3x) - 0.5\sin(5x), \quad
v^0 = 4\sin(x) - 1\sin(3x) +0.5\sin(5x).
$$

The same rule of indexsation applay to all vector coresponding to the finite part of points from the phase space.
- Parameters ${d_1,d_2,B,A}$ are taken from file `params.txt` in the form 

```
{0.2, 0.02, 2, 1}
-------------
{d_1,d_2,B,A}
```

##sampleDyn.cpp

Program for numerically integrating the Brusselator system. 

- Integrating in the space of odd coefficients.
- The initial data is taken from the content of file `initialValue.txt`. The initial data determines the size of Gallerkin projection.
- By using program:
- 
```
sampleDyn 100 10 1000
-------------
starting time, duration of simulation, number of steps
```
Will output uniformly sampled points $1000$ from the trajectory at times $[100,110]$ starting from initial date and save them to file `output2.txt`.

##makeSubsection.cpp

Program for finding periodic orbits of the Brusselator system. It tries to find a fixed periodic point by iterating the Galerkin approximation of some Poincare  map. Additionally after finding the point it proposed the section (additionally with intermediate sections if set).

- The initial data for searching is taken from the content of file `initialValue.txt`. This initial data determines the size of the Galerkin projection.
- In file `section.txt` we set the obtions for this program
```
0 0.7
1
1e-5
------
numVariable Value
numberOfSection
MaxAbsEigenValue
```
- $numVariable$ and $Value$ define the starting Poincare map and section for which we are searching for which we are searching for fixed point. The $numVariable$ defines the number of variable which is set to be equal to $Value$ at this Section.
- The $numberOfSection$ defines how many section program will output.
- The $MaxAbsEigenValue$ defines treshold for eigenvalues, for which the coresponding eigenvector are taken as the cordinates on ouputed section. Rest of them is taken by the other method.
- It the file 'subsections.txt' it ouputs subsections, which orgins are in trajectory of the numerical aproximation of periodic orbit.


## CAProof.cpp

This program is designed to prove the existence of a periodic orbit. It does so by checking whether for a defined set $X^0$ satisfies the following condition:

$$ P(X^0) \subset X^0 $$

For some Poincare map.
If this condition is satisfied, it confirms the existence of a periodic orbit. The program utilizes a rigorous C0 algorithm for integrating partial differential equations (PDEs) to compute the image.
Here are some details about the program and its inputs:

- The section which defined Poincare map is read from `subsections.txt`,
- In file `initialValuePoincare.txt`, we set:
```
{[-1e-4,1e-4],[-1e-4,1e-4],[-1e-4,1e-4],[-1e-7,1e-7]} {[-10,10],[4,4]}
-------------
L  {[C^-,C^+], [s,s]}

```
Which are used in defining $X^0$
- The set $X_0$ is defined as follows:
   $$X^0 = X_P^0 + X_Q^0,$$

where
$$X_P^0 = x^0+ Aq_0$$
 - Parameter $L$ defines entries of $q_0$: first is 0 the next are taken from the $L,$
and
$$
    X_Q^0 = (\sum_{i=n+1}^\infty \frac{[C^-,C^+]}{(2i-1)^s}\sin((2i-1)x),\sum_{i=n+1}^\infty \frac{[C^-,C^+]}{(2i-1)^s}\sin((2i-1)x) ),
$$
Here, $n$ depends on the size of the of size of setions in the 'subsection.txt', $C$ and $s$ from file `initialValuePoincare.txt`

   - The x^0 comes from orgin of Poincare section.
   - The number of modes explicitly represented.
   
The rigorous C0 integration is performed in the odd subspace of Fourier coefficients.

## Code Information

- The programs are using the [CAPD library](http://capd.ii.uj.edu.pl/index.php) - a tool for nonrigorous and validated numerics for dynamical systems.
  - Used version of the library: 5.1.2,
  - The Makefile assumes that the CAPD library is located in the following position relative to the main directory:
   ```
   # directory where capd scripts are (e.g. capd-config)
   CAPDBINDIR = ../capd-capdDynSys-5.1.2/bin/
   
   ```

- Folder DissipativePDE contains tools for rigorous integration of PDEs:
  - Folder DissipativePDE\Algebra contains the structure infinite series used in the algorithm and implementation of operations on them,
  - Folder DissipativePDE\Set contains structures for sets which are used in rigorous integration,
  - Folder DissipativePDE\VectorField contains structures for VectorFields used in rigorous integration,
  - Folder DissipativePDE\VectorFieldMaker contains additional methods that allow producing vector field fields for Gallerkin projection in string form, which is used in CAPD IMap class,
  - Folder DissipativePDE\SolverPDE contains methods and structs for rigorous C0 and C1 integrations.

- Folder Utils mainly contains Input/Output settings.

- Folder ChafeeInfante contains code specifically dedicated to Chafee-Infante:
   - Folder ChafeeInfante\ChafeeInfanteVecField containing C0 and C1 vector fields used in rigrous integration
   - Folder ChafeeInfante\GallerkinProjections containing method to produce strings, used in IMap class, for Gallerkin projection of the Chafee-Infante equation.





  



