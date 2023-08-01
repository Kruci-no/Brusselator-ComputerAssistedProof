# Brusselator-ComputerAssistedProof

The code is used for computer-assisted proof of the existence of a periodic orbit and its attraction for the non-autonomous Chafee-Infante system:

$$\begin{cases}
 u_t = d_1  u_{xx}- (B+1)u +u^2v + A \sin(x)   \text{ for } (x,t)\in (0,\pi)\times(0,\infty),\\
 v_t = d_2 v_{xx} + Bu - u^2v    \text{ for } (x,t)\in (0,\pi)\times(0,\infty),\\
u(t,x) = 0 \text{ for }  (x,t)\in \{0,\pi\} \times (0,\infty),\\
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
u^0 = 3\sin(x) - 0sin(3x) - 0.5\sin(5x),\\
v^0 = 4\sin(X) - 1\sin(3x) +0.5\sin(5x)
$$
The same rule of indexsation applay to all vector coresponding to the finite part of points from the phase space.
- Parameters ${d_1,d_2,B,A}$ are taken from file `params.txt` in the form 

```
{0.2, 0.02, 2, 1}
-------------
{d_1,d_2,B,A}
```

## ChafeeInfante/sampleDyn.cpp

Program for numerically integrating the Chafee-Infante system. 

- Integrating in the space of odd coefficients.
- The initial data is taken from the content of file `ChafeeInfante/textFiles/initialValue.txt`. The initial data determines the size of Gallerkin projection.
- File `ChafeeInfante/textFiles/sampleDynOptions.txt` contains options for the program:

```
0 1 2000
-------------
starting time, duration of simulation, number of steps
```

After running, the program should output sampled points from the trajectory and save them to file `ChafeeInfante/textFiles/sampleDynOutput.txt`.

## ChafeeInfante\findPeriodicPoint

Program for finding periodic orbits of the Chafee-Infate equation. It tries to find a fixed periodic point by iterating the Galerkin approximation of the map.

$$
T(u_0) = u(1,0;u_0).
$$

- The initial data for searching is taken from the content of file `ChafeeInfante\textFiles\initialValue.txt`. This initial data determines the size of the Galerkin projection.
- It assumes that parameter $\omega = 2\pi$, so the value of $\omega$ in file `ChafeeInfante\textFiles\params.txt` does not matter.
- It outputs the proposed approximation of the founded periodic orbit at zero. It also outputs an approximation of the variational matrix for the map $T$.


## ChafeeInfante\CAProof.cpp

This program is designed to prove the existence of a periodic orbit and demonstrate its local attraction. It does so by checking whether a defined set $X^0$ satisfies the following condition:

$$ T(X^0) \subset X^0 $$

If this condition is satisfied, it confirms the existence of a periodic orbit. The program utilizes a rigorous C0 algorithm for integrating partial differential equations (PDEs) to compute the image. Additionally, it attempts to prove that the orbit is locally attracting by verifying:

$$
||\frac{\partial T}{\partial x}(X_0)||_{C_0} < 1.
$$

The computation of derivatives employs a rigorous C1 integration algorithm.

Here are some details about the program and its inputs:

- The set $X^0$ is centered around an initial condition $u^0$ read from `ChafeeInfante\textFiles\initialValue.txt`,
- It assumes that parameter $\omega = 2\pi$, so the value of $\omega$ in file `ChafeeInfante\textFiles\params.txt` does not matter.
- Options for setting up the assisted proof can be found in `ChafeeInfante\textFiles\sampleDynOptions.txt`, including parameters such as:
```
1e-4 1 3 
6 14
9 17 0
-------------
eps , C, s
mainC0Size, fullC0Size
mainC1Size, fullC1Size, expColumns

```
- The set $X_0$ is defined as follows:
   $$X^0 = X_P^0 + X_Q^0,$$
  
where
$$X_P^0 = u^0(x)+ [-1,1]eps*\sum_{i=1}^n \sin((2i-1)x),$$

and

$$
    X_Q^0 = \sum_{i=n+1}^\infty \frac{C[-1,1]}{(2i-1)^s}\sin((2i-1)x).
$$

Here, $n$ depends on the size of the initial condition.

- For C0 computation, the parameters are set as follows:
   - `mainC0Size`: The number of modes considered for the differential inclusion.
   - `fullC0Size`: The number of modes explicitly represented.
   
The rigorous C0 integration is performed in the odd subspace of Fourier coefficients.

- For C1 computation, the parameters are set as follows:
   - `mainC1Size`: The number of modes considered for each variable $(u,h)$ in the differential inclusion.
   - `fullC1Size`: The number of modes explicitly represented for each variable $(u,h)$.
   - `expColumns`: Number of columns where derivatives are computed explicitly.

The rigorous C1 integration is conducted in full space, not just in the odd subspace. In this case, there are variables representing variational equations. Therefore, the actual number of modes used for differential inclusion is $2 * \text{mainC1Size}$ and the number of modes explicitly represented is $2 * \text{fullC1Size}$.

The remaining columns are computed by setting them to be:

$$
    X^{0,h} = \sum_{i=\text{expColumns}+1}^\infty u_i\sin(ix),\quad \text{where } u_i\in[-1,1].
$$

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





  



