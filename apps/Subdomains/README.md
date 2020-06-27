# Frac3D/Subdomains

Solution of 3D fluid flow problems in porous media with three subdomains.

## Problem Description

The problem is to find <a href="https://www.codecogs.com/eqnedit.php?latex=p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p" title="p" /></a> such that

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;-\nabla\cdot\left(K\nabla&space;p\right)=0,&space;&\;\text{in}\;\Omega,&space;\\&space;p=\overline{p},&space;&\;\text{in}\;\partial\Omega_p,&space;\\&space;-K\nabla&space;p\cdot&space;n=0,&space;&\;\text{in}\;\partial\Omega_u,&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;-\nabla\cdot\left(K\nabla&space;p\right)=0,&space;&\;\text{in}\;\Omega,&space;\\&space;p=\overline{p},&space;&\;\text{in}\;\partial\Omega_p,&space;\\&space;-K\nabla&space;p\cdot&space;n=0,&space;&\;\text{in}\;\partial\Omega_u,&space;\end{align*}" title="\begin{align*} -\nabla\cdot\left(K\nabla p\right)=0, &\;\text{in}\;\Omega, \\ p=\overline{p}, &\;\text{in}\;\partial\Omega_p, \\ -K\nabla p\cdot n=0, &\;\text{in}\;\partial\Omega_u, \end{align*}" /></a>

and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\textstyle\partial\Omega_p\cap\partial\Omega_u=\emptyset" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\textstyle\partial\Omega_p\cap\partial\Omega_u=\emptyset" title="\textstyle\partial\Omega_p\cap\partial\Omega_u=\emptyset" /></a>.

The weak formulation is find <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\textstyle&space;p\in&space;V(\Omega)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\textstyle&space;p\in&space;V(\Omega)" title="\textstyle p\in V(\Omega)" /></a>, such that

<a href="https://www.codecogs.com/eqnedit.php?latex=\int_{\Omega}K\nabla&space;p\cdot\nabla&space;v\;\text{d}x=0," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\int_{\Omega}K\nabla&space;p\cdot\nabla&space;v\;\text{d}x=0," title="\int_{\Omega}K\nabla p\cdot\nabla v\;\text{d}x=0," /></a>

for all <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;v\in&space;V(\Omega)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;v\in&space;V(\Omega)" title="v\in V(\Omega)" /></a>.

## Numerical Solution

### Subdomains:

The domain is decomposed into three subdomains:
- Upper subdomain <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Omega_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Omega_1" title="\Omega_1" /></a> (gray color, ID 19), with <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;K_1=10^{-6}\;\text{m}^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;K_1=10^{-6}\;\text{m}^2" title="K_1=10^{-6}\;\text{m}^2" /></a>;
- Lower subdomain <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Omega_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Omega_2" title="\Omega_2" /></a> (red color, ID 20), with <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;K_2=10^{-6}\;\text{m}^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;K_2=10^{-6}\;\text{m}^2" title="K_2=10^{-6}\;\text{m}^2" /></a>;
- Bottom subdomain <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Omega_3" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Omega_3" title="\Omega_3" /></a> (blue color, ID 18), with <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;K_3=10^{-5}\;\text{m}^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;K_3=10^{-5}\;\text{m}^2" title="K_3=10^{-5}\;\text{m}^2" /></a>.

<p float="left">
	<img src="figs/Subdomains.png" alt="subdomains" height=300/>
</p>

### Boundaries:

The boundaries are subdivided into two Dirichlet boundaries and one Neumann boundary:
- Dirichlet boundary <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial\Omega_{p1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\partial\Omega_{p1}" title="\partial\Omega_{p1}" /></a> (narrow band defined by <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\left(0\;\text{m},&space;0\;\text{m}\right)\times\left(0\;\text{m},&space;100\;\text{m}\right)\times\left(90\;\text{m},&space;100\;\text{m}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\left(0\;\text{m},&space;0\;\text{m}\right)\times\left(0\;\text{m},&space;100\;\text{m}\right)\times\left(90\;\text{m},&space;100\;\text{m}\right)" title="\left(0\;\text{m}, 0\;\text{m}\right)\times\left(0\;\text{m}, 100\;\text{m}\right)\times\left(90\;\text{m}, 100\;\text{m}\right)" /></a>, ID 32), with prescribed <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;p=4\;\text{m}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;p=4\;\text{m}" title="p=4\;\text{m}" /></a>;
- Dirichlet boundary <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial\Omega_{p2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\partial\Omega_{p2}" title="\partial\Omega_{p2}" /></a> (narrow band defined by <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\left(0\;\text{m},&space;100\;\text{m}\right)\times\left(0\;\text{m},&space;0\;\text{m}\right)\times\left(0\;\text{m},&space;10\;\text{m}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\left(0\;\text{m},&space;100\;\text{m}\right)\times\left(0\;\text{m},&space;0\;\text{m}\right)\times\left(0\;\text{m},&space;10\;\text{m}\right)" title="\left(0\;\text{m}, 100\;\text{m}\right)\times\left(0\;\text{m}, 0\;\text{m}\right)\times\left(0\;\text{m}, 10\;\text{m}\right)" /></a>, ID 31), with prescribed <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;p=1\;\text{m}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;p=1\;\text{m}" title="p=1\;\text{m}" /></a>;
- Neumann boundary <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\partial\Omega_u=\partial\Omega-\left(\partial\Omega_{p1}\cup\partial\Omega_{p2}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\partial\Omega_u=\partial\Omega-\left(\partial\Omega_{p1}\cup\partial\Omega_{p2}\right)" title="\partial\Omega_u=\partial\Omega-\left(\partial\Omega_{p1}\cup\partial\Omega_{p2}\right)" /></a>.

<p float="left">
	<img src="figs/Boundaries.png" alt="boundaries" height=300/>
</p>

### Results:

The solution is obtained using second degree CG elements for the pressure field and zero-order DG elements for the permeability field.

<p float="left">
	<img src="figs/Solution.png" alt="solution" height=300/>
</p>