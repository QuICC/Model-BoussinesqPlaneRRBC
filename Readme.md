# Boussinesq Rotating Rayleigh-Benard convection(RRBC) in a plane layer

The model equations are

$$
(\partial_t -\Delta){\bf u}= {\bf u}\times(\nabla\times{\bf u})- \frac{1}{E}\hat{\bf z}\times{\bf u}  -\nabla \Pi + \frac{Ra}{Pr} \Theta\hat{\bf z}\\
$$

$$
(\partial_t -\frac{1}{Pr}\Delta)\Theta = -{\bf u}\cdot\nabla\Theta 
$$

with the parameters defined as

$$
E = \frac{\nu}{2 \Omega H^2}
$$

$$
Ra = \frac{\alpha g \beta H^4}{\nu \kappa}
$$

$$
Pr = \frac{\nu}{\kappa}
$$

where $H$ is the height of the layer (considered doubly perioxdic along the horizontal directions) and $\beta$ is the super-adiabatic temperature gradient.

To solve the system using a Toroidal/Poloidal decomposition of the velocity

$$
{\bf u}(x,y,z) = \nabla\times T(x,y,z) \hat{\bf z} + \nabla\times\nabla\times P(x,y,z) \hat{\bf z} + \overline{u_x}(z) \hat{\bf x} + \overline{u_y}(z) \hat{\bf y}
$$

where the overline denote horizontally averaged quantities that depend only on the vertical coordinate, $z$. A similar decomposition is made for temperature:

$$
\Theta(x,y,z) = \overline{\theta}(z) + \theta(x,y,z)
$$