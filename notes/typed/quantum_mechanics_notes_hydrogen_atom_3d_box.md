# Introductory Quantum Mechanics: Derivations & Notes

> **Organization Note:** These notes are transcribed directly from your original pages, preserving your exact variable choices (such as $y(r) = rR(r)$, $x = r\epsilon$, $j = n - l - 1$, and $D = \frac{d}{dx}$), intermediate steps, and worked examples. The sections are ordered logically from the 3D Box to the Angular Momentum Operator, followed by the full Hydrogen Atom derivation.

## How these notes connect to the code

These derivations are the theory context for the project, not a description of a single end-to-end implementation:

* **Part 1, 3-D box:** there is no matching 3-D box solver in this workspace. [Finite Potential Well Solutions.py](../../src/Finite%20Potential%20Well%20Solutions.py) explores a separate one-dimensional finite well.
* **Part 2, angular momentum:** this is a derivation only; the project does not implement a general angular-momentum operator.
* **Part 3, hydrogen atom:** [Hydrogen Electron Probability Cloud.py](../../src/Hydrogen%20Electron%20Probability%20Cloud.py) plots selected states using hand-written expressions. It does not calculate the general radial function and spherical harmonics derived below, and its visual output has not been validated against these notes.

For known implementation caveats, see the [errata](../../ERRATA.md) and [archive README](../../README.md).

---

## Part 1: Particle in a 3-D Box *(Original Pages 23–29)*

### 1.1 System Setup & Potential
Consider a particle inside a 3-D cubic box of dimensions $L_x = L_y = L_z = L$:
* **Potential inside the box:** $V(x,y,z) = 0$ when $0 < x < L$, $0 < y < L$, and $0 < z < L$.
* **Potential outside the box:** Infinite ($\infty$), meaning the wave function at the walls is zero.

The Time-Independent Schrödinger Equation (TISE) in 1-D is:
$$-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2} + V(x)\psi(x) = E\psi(x)$$

In 3-D, the wave function depends on three spatial coordinates, $\psi(x,y,z)$, and since $V(x,y,z) = 0$ inside the box, the 3-D Schrödinger equation becomes:
$$-\frac{\hbar^2}{2m}\left(\frac{\partial^2\psi(x,y,z)}{\partial x^2} + \frac{\partial^2\psi(x,y,z)}{\partial y^2} + \frac{\partial^2\psi(x,y,z)}{\partial z^2}\right) = E\psi(x,y,z)$$

### 1.2 Separation of Variables
Assume a separable solution:
$$\psi(x,y,z) = X(x)Y(y)Z(z)$$

This separates the PDE into three ordinary differential equations (ODEs), where the total energy is the sum of the energies of the individual components, $E = \epsilon_x + \epsilon_y + \epsilon_z$:
$$-\frac{\hbar^2}{2m}\frac{d^2X(x)}{dx^2} = \epsilon_x X(x), \qquad -\frac{\hbar^2}{2m}\frac{d^2Y(y)}{dy^2} = \epsilon_y Y(y), \qquad -\frac{\hbar^2}{2m}\frac{d^2Z(z)}{dz^2} = \epsilon_z Z(z)$$

### 1.3 Solving the ODEs
For the $x$-component:
$$\frac{d^2X(x)}{dx^2} = -\frac{2m\epsilon_x}{\hbar^2}X(x)$$

Let $k_x^2 = \frac{2m\epsilon_x}{\hbar^2}$:
$$\frac{d^2X(x)}{dx^2} + k_x^2 X(x) = 0$$

Using the differential operator $D = \frac{d}{dx}$, we write the auxiliary equation ($A.E.$):
$$(D^2 + k_x^2)X(x) = 0 \implies m^2 + k_x^2 = 0 \implies m = \pm ik_x$$

Since the roots are complex, the general solution takes the form $y = e^{\alpha x}[c_1\cos(\beta x) + c_2\sin(\beta x)]$ with $\alpha = 0$:
$$X(x) = C_1\cos(k_x x) + C_2\sin(k_x x)$$

Relabeling constants $C_2 = A_x$ and $C_1 = B_x$ for convenience, the solutions for all three components are:
$$X(x) = A_x\sin(k_x x) + B_x\cos(k_x x)$$
$$Y(y) = A_y\sin(k_y y) + B_y\cos(k_y y)$$
$$Z(z) = A_z\sin(k_z z) + B_z\cos(k_z z)$$

Expanding $\psi(x,y,z) = X(x)Y(y)Z(z)$:
$$\psi(x,y,z) = \big(A_x\sin(k_x x) + B_x\cos(k_x x)\big)\big(A_y\sin(k_y y) + B_y\cos(k_y y)\big)\big(A_z\sin(k_z z) + B_z\cos(k_z z)\big)$$

### 1.4 Applying Boundary Conditions & Quantized Energy
1. **Boundary condition at origin ($x=0, y=0, z=0$):**
   $$\psi(0,y,z) = 0, \qquad \psi(x,0,z) = 0, \qquad \psi(x,y,0) = 0$$
   Substituting these in yields $B_x = B_y = B_z = 0$. Letting $A = A_x A_y A_z$ and expanding $k_x, k_y, k_z$:
   $$\psi(x,y,z) = A\sin\left(\sqrt{\frac{2m\epsilon_x}{\hbar^2}}\,x\right)\sin\left(\sqrt{\frac{2m\epsilon_y}{\hbar^2}}\,y\right)\sin\left(\sqrt{\frac{2m\epsilon_z}{\hbar^2}}\,z\right)$$

2. **Boundary condition at distance $L$ ($x=L, y=L, z=L$):**
   $$\psi(L,y,z) = 0, \qquad \psi(x,L,z) = 0, \qquad \psi(x,y,L) = 0$$
   Since the other two spatial coordinates may be non-zero in each case, we cannot depend on their sine terms to vanish the equation. Thus:
   $$\sin\left(\sqrt{\frac{2m\epsilon_x}{\hbar^2}}\,L\right) = 0, \qquad \sin\left(\sqrt{\frac{2m\epsilon_y}{\hbar^2}}\,L\right) = 0, \qquad \sin\left(\sqrt{\frac{2m\epsilon_z}{\hbar^2}}\,L\right) = 0$$

   Solving the $x$-component (noting $n_x \neq 0$, otherwise $\psi = 0$):
   $$\sqrt{\frac{2m\epsilon_x}{\hbar^2}}\,L = n_x\pi \implies \frac{2m\epsilon_x L^2}{\hbar^2} = n_x^2\pi^2 \implies \epsilon_x = \frac{n_x^2\pi^2\hbar^2}{2mL^2} \quad (n_x = 1, 2, 3, \dots)$$

   Similarly, for $y$ and $z$:
   $$\epsilon_y = \frac{n_y^2\pi^2\hbar^2}{2mL^2}, \qquad \epsilon_z = \frac{n_z^2\pi^2\hbar^2}{2mL^2} \quad (n_y, n_z = 1, 2, 3, \dots)$$

   **Total Energy:**
   $$E = \epsilon_x + \epsilon_y + \epsilon_z = \frac{\hbar^2\pi^2}{2mL^2}\left(n_x^2 + n_y^2 + n_z^2\right)$$

   Substituting $\epsilon_x, \epsilon_y, \epsilon_z$ back into $\psi(x,y,z)$ simplifies the wave function to:
   $$\psi(x,y,z) = A\sin\left(\frac{n_x\pi x}{L}\right)\sin\left(\frac{n_y\pi y}{L}\right)\sin\left(\frac{n_z\pi z}{L}\right)$$

### 1.5 Normalization
Using the normalization condition $\int_0^L \int_0^L \int_0^L \psi\psi^*\,dx\,dy\,dz = 1$:
$$A^2 \int_0^L \sin^2\left(\frac{n_x\pi x}{L}\right)dx \int_0^L \sin^2\left(\frac{n_y\pi y}{L}\right)dy \int_0^L \sin^2\left(\frac{n_z\pi z}{L}\right)dz = 1$$

Solving the $x$-integral via substitution $u_x = \frac{n_x\pi x}{L}$ ($dx = \frac{L\,du_x}{n_x\pi}$, limits $0 \to n_x\pi$) and using $\sin^2 u = \frac{1 - \cos 2u}{2}$:
$$\int_0^L \sin^2\left(\frac{n_x\pi x}{L}\right)dx = \frac{L}{n_x\pi}\int_0^{n_x\pi} \frac{1 - \cos(2u_x)}{2}\,du_x = \frac{L}{2n_x\pi}\left[n_x\pi - \int_0^{n_x\pi}\cos(2u_x)\,du_x\right]$$

Substituting $v_x = 2u_x$ ($du_x = \frac{dv_x}{2}$, limits $0 \to 2n_x\pi$):
$$= \frac{L}{2n_x\pi}\left[n_x\pi - \frac{1}{2}\Big[\sin v_x\Big]_0^{2n_x\pi}\right] = \frac{L}{2n_x\pi}\left[n_x\pi - \frac{1}{2}(\sin(2n_x\pi) - \sin(0))\right] = \frac{L}{2n_x\pi}(n_x\pi) = \frac{L}{2}$$

Solving identically for the $y$ and $z$ integrals gives:
$$A^2 \left(\frac{L}{2}\right)\left(\frac{L}{2}\right)\left(\frac{L}{2}\right) = 1 \implies \frac{A^2 L^3}{8} = 1 \implies A = \sqrt{\frac{8}{L^3}}$$

**Normalized 3-D Wave Function:**
* For a cubic box ($L_x = L_y = L_z = L$, volume $V = L^3$):
  $$\psi(x,y,z) = \sqrt{\frac{8}{V}}\sin\left(\frac{n_x\pi x}{L}\right)\sin\left(\frac{n_y\pi y}{L}\right)\sin\left(\frac{n_z\pi z}{L}\right)$$
* For a general cuboid box ($L_x \neq L_y \neq L_z$, volume $V = L_x L_y L_z$):
  $$\psi(x,y,z) = \sqrt{\frac{8}{V}}\sin\left(\frac{n_x\pi x}{L_x}\right)\sin\left(\frac{n_y\pi y}{L_y}\right)\sin\left(\frac{n_z\pi z}{L_z}\right)$$
  where $n_x, n_y, n_z \in \mathbb{N}$.

---

## Part 2: Angular Momentum Operator in Spherical Coordinates *(Original Pages 19–21)*

### 2.1 Classical to Quantum Cartesian Components
In classical mechanics, orbital (revolution) angular momentum for a point mass $m$ is the cross product $\vec{L} = \vec{r} \times \vec{p}$:
$$L_x = y p_z - z p_y, \qquad L_y = z p_x - x p_z, \qquad L_z = x p_y - y p_x$$

Replacing momentum components with the quantum mechanical momentum operator $\hat{p} = -i\hbar\nabla = -i\hbar\left(\hat{i}\frac{\partial}{\partial x} + \hat{j}\frac{\partial}{\partial y} + \hat{k}\frac{\partial}{\partial z}\right)$:
$$\hat{L}_x = -i\hbar\left(y\frac{\partial}{\partial z} - z\frac{\partial}{\partial y}\right)$$
$$\hat{L}_y = -i\hbar\left(z\frac{\partial}{\partial x} - x\frac{\partial}{\partial z}\right)$$
$$\hat{L}_z = -i\hbar\left(x\frac{\partial}{\partial y} - y\frac{\partial}{\partial x}\right)$$

### 2.2 Conversion to Spherical Coordinates & Deriving $\hat{L}^2$
Converting the Cartesian operators into spherical coordinates $(r, \theta, \phi)$:
$$\hat{L}_x = i\hbar\left(\sin\phi\frac{\partial}{\partial\theta} + \cot\theta\cos\phi\frac{\partial}{\partial\phi}\right)$$
$$\hat{L}_y = i\hbar\left(-\cos\phi\frac{\partial}{\partial\theta} + \cot\theta\sin\phi\frac{\partial}{\partial\phi}\right)$$
$$\hat{L}_z = -i\hbar\frac{\partial}{\partial\phi}$$

To find $\hat{L}^2 = \hat{L}_x^2 + \hat{L}_y^2 + \hat{L}_z^2$, we square each component (applying the product rule carefully to the mixed angular derivatives):
$$\hat{L}_x^2 = -\hbar^2\left(\sin^2\phi\frac{\partial^2}{\partial\theta^2} + \cot^2\theta\cos^2\phi\frac{\partial^2}{\partial\phi^2} + \sin\phi\cos\phi\cot\theta\frac{\partial}{\partial\theta}\frac{\partial}{\partial\phi} + \dots\right)$$
$$\hat{L}_y^2 = -\hbar^2\left(\cos^2\phi\frac{\partial^2}{\partial\theta^2} + \cot^2\theta\sin^2\phi\frac{\partial^2}{\partial\phi^2} - \cos\phi\sin\phi\cot\theta\frac{\partial}{\partial\theta}\frac{\partial}{\partial\phi} + \dots\right)$$
$$\hat{L}_z^2 = -\hbar^2\frac{\partial^2}{\partial\phi^2}$$

Adding all three terms together (using $\sin^2\phi + \cos^2\phi = 1$ and $\cot^2\theta + 1 = \csc^2\theta = \frac{1}{\sin^2\theta}$):
$$\hat{L}^2 = -\hbar^2\left[\frac{\partial^2}{\partial\theta^2} + \cot\theta\frac{\partial}{\partial\theta} + \left(\cot^2\theta + 1\right)\frac{\partial^2}{\partial\phi^2}\right]$$
$$\hat{L}^2 = -\hbar^2\left[\frac{\partial^2}{\partial\theta^2} + \frac{\cos\theta}{\sin\theta}\frac{\partial}{\partial\theta} + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right]$$

Combining the $\theta$-derivatives via the product rule $\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) = \frac{\partial^2}{\partial\theta^2} + \cot\theta\frac{\partial}{\partial\theta}$, we obtain:
$$\hat{L}^2 = -\hbar^2\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right]$$

---

## Part 3: Quantum Mechanics of the Hydrogen Atom *(Original Pages 1, 3–18)*

### 3.1 Physical Setup & Hamiltonian in Spherical Coordinates
Consider an electron (charge $q_1 = -e$, mass $m_e$) orbiting a proton nucleus (charge $q_2 = +e$, mass $m_p$) at position vector $r$, inclination (polar) angle $\theta$, and azimuth angle $\phi$.
* **Reduced mass:** $\mu = \frac{m_e m_p}{m_e + m_p}$
* **Coulomb potential:**
  $$V(r) = \frac{k q_1 q_2}{r} = -\frac{e^2}{4\pi\epsilon_0 r} \quad \text{where } k = \frac{1}{4\pi\epsilon_0}$$

Starting from the 3-D Time-Independent Schrödinger Equation in Cartesian coordinates:
$$-\frac{\hbar^2}{2\mu}\nabla^2\psi(x,y,z) + V(x,y,z)\psi(x,y,z) = E\psi(x,y,z)$$

Using the spherical coordinate transformations $x = r\sin\theta\cos\phi$, $y = r\sin\theta\sin\phi$, and $z = r\cos\theta$, the Laplacian operator in spherical coordinates is:
$$\nabla^2 = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2}{\partial\phi^2}$$

Thus, the Hamiltonian operator $\hat{H}(r,\theta,\phi)$ becomes:
$$\hat{H}(r,\theta,\phi) = -\frac{\hbar^2}{2\mu r^2}\left[\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right) + \frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right] - \frac{e^2}{4\pi\epsilon_0 r}$$

### 3.2 Separation of Radial and Angular Variables
Writing the Schrödinger equation $\hat{H}\psi(r,\theta,\phi) = E\psi(r,\theta,\phi)$ and multiplying both sides by $2\mu r^2$:
$$-\hbar^2\frac{\partial}{\partial r}\left(r^2\frac{\partial\psi}{\partial r}\right) - \hbar^2\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial\psi}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2\psi}{\partial\phi^2}\right] - \frac{2\mu r^2 e^2}{4\pi\epsilon_0 r}\psi = 2\mu r^2 E\psi$$

Recognizing the angular momentum squared operator $\hat{L}^2$ in the second term and moving it to the right-hand side:
$$\underbrace{\hbar^2\frac{\partial}{\partial r}\left(r^2\frac{\partial\psi(r,\theta,\phi)}{\partial r}\right) + 2\mu r^2\big(E - V(r)\big)\psi(r,\theta,\phi)}_{\text{Radial motion}} = \underbrace{\hat{L}^2\psi(r,\theta,\phi)}_{\text{Rotational motion}}$$

Because $\hat{L}^2$ does not act on $r$, we apply separation of variables with $\psi(r,\theta,\phi) = R(r)Y(\theta,\phi)$ (where $R(r)$ gives the radial distance from the proton and $Y(\theta,\phi)$ is the spherical harmonic function giving angular position):
$$\hbar^2 Y(\theta,\phi)\frac{d}{dr}\left(r^2\frac{dR(r)}{dr}\right) + 2\mu r^2\big(E - V(r)\big)R(r)Y(\theta,\phi) = R(r)\hat{L}^2 Y(\theta,\phi)$$

Dividing both sides by $R(r)Y(\theta,\phi)$ and setting both sides equal to the separation constant $\lambda = l(l+1)\hbar^2$ (where $l$ is the **azimuthal/orbital angular momentum quantum number**):

1. **Radial Equation (1):**
   $$\frac{\hbar^2}{R(r)}\frac{d}{dr}\left(r^2\frac{dR(r)}{dr}\right) + 2\mu r^2\big(E - V(r)\big) = l(l+1)\hbar^2$$
   Dividing by $\hbar^2$:
   $$\frac{1}{R(r)}\frac{d}{dr}\left(r^2\frac{dR(r)}{dr}\right) - \frac{2\mu r^2}{\hbar^2}\big(V(r) - E\big) = l(l+1)$$

2. **Angular Equation (2):**
   $$\frac{1}{Y(\theta,\phi)}\hat{L}^2 Y(\theta,\phi) = l(l+1)\hbar^2$$
   Expanding $\hat{L}^2$ and dividing by $-\hbar^2$:
   $$\frac{1}{Y(\theta,\phi)\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial Y(\theta,\phi)}{\partial\theta}\right) + \frac{1}{Y(\theta,\phi)\sin^2\theta}\frac{\partial^2 Y(\theta,\phi)}{\partial\phi^2} = -l(l+1)$$

---

### 3.3 Solving the Angular Equation
Separate the spherical harmonic function as $Y(\theta,\phi) = f(\theta)g(\phi)$. Substituting into the angular equation and multiplying both sides by $\sin^2\theta$:
$$\frac{\sin\theta}{f(\theta)}\frac{d}{d\theta}\left(\sin\theta\frac{df(\theta)}{d\theta}\right) + l(l+1)\sin^2\theta = -\frac{1}{g(\phi)}\frac{d^2 g(\phi)}{d\phi^2}$$

Equating both sides to the separation constant $m^2$ (where $m$ is the **magnetic quantum number**, describing orbital orientation):
$$\frac{\sin\theta}{f(\theta)}\frac{d}{d\theta}\left(\sin\theta\frac{df(\theta)}{d\theta}\right) + l(l+1)\sin^2\theta = m^2 \qquad \text{(Equation 3: Polar)}$$
$$\frac{1}{g(\phi)}\frac{d^2 g(\phi)}{d\phi^2} = -m^2 \qquad \text{(Equation 4: Azimuthal)}$$

#### A. Azimuthal Angle Solution ($g(\phi)$)
$$\frac{d^2 g(\phi)}{d\phi^2} + m^2 g(\phi) = 0 \implies \left(\frac{d}{d\phi} + im\right)\left(\frac{d}{d\phi} - im\right)g(\phi) = 0$$

Taking the branch $\frac{dg(\phi)}{d\phi} - im g(\phi) = 0$:
$$\frac{dg(\phi)}{g(\phi)} = im\,d\phi \implies \int \frac{1}{g(\phi)}\,dg(\phi) = \int im\,d\phi \implies \ln\{g(\phi)\} = im\phi$$
$$g_m(\phi) = e^{im\phi}$$

#### B. Polar Angle Solution ($f(\theta)$)
Multiplying Equation (3) through by $f(\theta)$ and expanding the derivative:
$$\sin^2\theta\frac{d^2 f(\theta)}{d\theta^2} + \sin\theta\cos\theta\frac{df(\theta)}{d\theta} + l(l+1)\sin^2\theta\,f(\theta) - m^2 f(\theta) = 0$$

Let $x = \cos\theta$, so $\frac{dx}{d\theta} = -\sin\theta$. Using the chain rule:
$$\frac{df(\theta)}{d\theta} = -\sin\theta\frac{df(x)}{dx}$$
$$\frac{d^2 f(\theta)}{d\theta^2} = \frac{d}{d\theta}\left(-\sin\theta\frac{df(x)}{dx}\right) = \sin^2\theta\frac{d^2 f(x)}{dx^2} - \cos\theta\frac{df(x)}{dx}$$

Substituting these derivatives back into the polar equation and dividing by $\sin^2\theta$ (with $\sin^2\theta = 1 - \cos^2\theta = 1 - x^2$):
$$(1 - x^2)\frac{d^2 f(x)}{dx^2} - 2x\frac{df(x)}{dx} + \left[l(l+1) - \frac{m^2}{1 - x^2}\right]f(x) = 0$$

This is the **Associated Legendre Differential Equation**. Its solutions are the **Associated Legendre Polynomials**:
$$P_l^m(x) = (-1)^m (1 - x^2)^{m/2}\frac{d^m}{dx^m}\big(P_l(x)\big)$$
where $P_l(x)$ are the ordinary Legendre polynomials given by Rodrigues' formula:
$$P_l(x) = \frac{(-1)^l}{2^l l!}\frac{d^l}{dx^l}(1 - x^2)^l$$

#### C. Normalized Spherical Harmonics ($Y_l^m(\theta,\phi)$)
Combining $f(\theta)$ and $g(\phi)$ with the normalization constant gives the orthonormal spherical harmonics (for $m \ge 0$):
$$Y_l^m(\theta,\phi) = (-1)^m \sqrt{\frac{(2l+1)(l-m)!}{4\pi(l+m)!}}\,P_l^m(\cos\theta)e^{im\phi}$$
*(For $m < 0$, $Y_l^{-|m|}(\theta,\phi) = (-1)^m \big(Y_l^{|m|}(\theta,\phi)\big)^*$.)*

---

### 3.4 Solving the Radial Equation & Quantized Energy
Expressing the potential in CGS units as $V(r) = -\frac{e^2}{r}$ (converting back to SI via $e^2 \to \frac{e^2}{4\pi\epsilon_0}$ at the end) and multiplying the radial equation by $R(r)$:
$$\frac{d}{dr}\left(r^2\frac{dR(r)}{dr}\right) + \left[\frac{2\mu r^2 e^2}{\hbar^2 r} + \frac{2\mu r^2 E}{\hbar^2} - l(l+1)\right]R(r) = 0$$

Let $y(r) = rR(r)$, so $R(r) = \frac{y(r)}{r}$. Evaluating the radial derivative term:
$$\frac{d}{dr}\left(r^2\frac{d}{dr}(r^{-1}y(r))\right) = \frac{d}{dr}\Big(r^2\big[-r^{-2}y(r) + r^{-1}y'(r)\big]\Big) = \frac{d}{dr}\big[-y(r) + ry'(r)\big] = r\frac{d^2 y(r)}{dr^2}$$

Substituting $R(r) = \frac{y(r)}{r}$ into the radial equation and dividing through by $r$:
$$\frac{d^2 y(r)}{dr^2} + \left[\frac{2\mu e^2}{\hbar^2 r} + \frac{2\mu E}{\hbar^2} - \frac{l(l+1)}{r^2}\right]y(r) = 0$$

Introduce the substitutions $\left(\frac{\epsilon}{2}\right)^2 = -\frac{2\mu E}{\hbar^2}$ and $x = r\epsilon$ (so $r = \frac{x}{\epsilon}$ and $\frac{d^2 y}{dr^2} = \epsilon^2\frac{d^2 y}{dx^2}$), then divide through by $\epsilon^2$:
$$\frac{d^2 y(x)}{dx^2} + \left[-\frac{1}{4} + \frac{2\mu e^2}{\hbar^2\epsilon x} - \frac{l(l+1)}{x^2}\right]y(x) = 0$$

Comparing this to the **Associated Laguerre Differential Equation**:
$$y_j^k(x)'' + \left[-\frac{1}{4} + \frac{2j + k + 1}{2x} - \frac{k^2 - 1}{4x^2}\right]y_j^k(x) = 0$$
whose solution is:
$$y_j^k(x) = e^{-x/2} x^{(k+1)/2} L_j^k(x)$$
where the Associated Laguerre Polynomial $L_j^k(x)$ and Laguerre Polynomial $L_j(x)$ are:
$$L_j^k(x) = (-1)^k \frac{d^k}{dx^k}L_{j+k}(x), \qquad L_j(x) = e^x \frac{d^j}{dx^j}\left(e^{-x}x^j\right)$$

Matching coefficients between the two equations:
1. $l(l+1) = \frac{k^2 - 1}{4} \implies k^2 = 4l^2 + 4l + 1 = (2l+1)^2 \implies k = 2l + 1$
2. $\frac{2\mu e^2}{\hbar^2\epsilon} = \frac{2j + k + 1}{2} = \frac{2j + 2l + 2}{2} = j + l + 1$

Since the indices of Laguerre polynomials must be non-negative integers, $j + l + 1$ must be a natural number. Define the **Principal Quantum Number** $n = j + l + 1$ (so $j = n - l - 1$).

#### Deriving Bohr Radius ($a_0$) & Energy ($E_n$)
Using the Bohr radius $a_0 = \frac{\hbar^2}{\mu e^2} = 0.529\text{ \AA} = 0.0529\text{ nm}$ (in SI: $a_0 = \frac{4\pi\epsilon_0\hbar^2}{\mu e^2}$) and $\frac{2\mu e^2}{\hbar^2\epsilon} = n$:
$$\left(\frac{2}{\epsilon}\right)^2 \left(\frac{\mu e^2}{\hbar^2}\right)^2 = n^2 \implies \left(-\frac{\hbar^2}{2\mu E}\right)\left(\frac{1}{a_0}\right)^2 = n^2$$

Solving for the quantized energy $E_n$:
$$E_n = -\frac{\hbar^2}{2\mu n^2 a_0^2}$$

#### Complete Radial & Total Wave Function
Substituting $k = 2l + 1$, $j = n - l - 1$, and $x = r\epsilon = \frac{2r}{a_0 n}$ into $y_j^k(x)$:
$$y_n^l(r) = e^{-r/(na_0)}\left(\frac{2r}{a_0 n}\right)^{l+1} L_{n-l-1}^{2l+1}\left(\frac{2r}{a_0 n}\right)$$

Since $R_n^l(r) = \frac{y_n^l(r)}{r}$, factoring out $r$ and including the radial normalization constant yields:
$$R_n^l(r) = \sqrt{\left(\frac{2}{a_0 n}\right)^3 \frac{(n - l - 1)!}{2n\big[(n + l)!\big]^3}}\,e^{-r/(na_0)}\left(\frac{2r}{a_0 n}\right)^l L_{n-l-1}^{2l+1}\left(\frac{2r}{a_0 n}\right)$$

**Complete Hydrogen Atom Wave Function ($\psi_{n,l,m} = R_n^l(r)Y_l^m(\theta,\phi)$):**
$$\psi_{n,l,m}(r,\theta,\phi) = \sqrt{\left(\frac{2}{a_0 n}\right)^3 \frac{(n - l - 1)!}{2n\big[(n + l)!\big]^3}}\,e^{-r/(na_0)}\left(\frac{2r}{a_0 n}\right)^l L_{n-l-1}^{2l+1}\left(\frac{2r}{a_0 n}\right) \cdot (-1)^m\sqrt{\frac{(2l+1)(l-m)!}{4\pi(l+m)!}}\,P_l^m(\cos\theta)e^{im\phi}$$

---

### 3.5 Worked Example: Computing $\psi_{211}(r,\theta,\phi)$ *(Original Pages 17–18)*
For the state $(n, l, m) = (2, 1, 1)$:

#### Step 1: Radial Part $R_2^1(r)$
$$R_2^1(r) = \sqrt{\left(\frac{2}{2a_0}\right)^3 \frac{(2 - 1 - 1)!}{2(2)\big[(2 + 1)!\big]^3}}\,e^{-r/(2a_0)}\left(\frac{2r}{2a_0}\right)^1 L_{2-1-1}^{2(1)+1}\left(\frac{2r}{2a_0}\right)$$
$$R_2^1(r) = \sqrt{\frac{1}{a_0^3} \cdot \frac{0!}{4(3!)^3}}\,e^{-r/(2a_0)}\left(\frac{r}{a_0}\right) L_0^3\left(\frac{r}{a_0}\right)$$

1. **Normalization prefactor:**
   $$\sqrt{\frac{1}{a_0^3 \cdot 4 \cdot 6^3}} = \sqrt{\frac{1}{864\,a_0^3}} = \frac{1}{12\sqrt{6}\,a_0^{3/2}}$$

2. **Associated Laguerre Polynomial $L_0^3(x)$:**
   First compute $L_3(x)$:
   $$L_3(x) = e^x \frac{d^3}{dx^3}\left(e^{-x}x^3\right) = e^x e^{-x}\left(-x^3 + 9x^2 - 18x + 6\right) = -x^3 + 9x^2 - 18x + 6$$
   Then differentiate 3 times for $L_0^3(x)$:
   $$L_0^3(x) = (-1)^3 \frac{d^3}{dx^3}L_3(x) = -\frac{d^2}{dx^2}\left(-3x^2 + 18x - 18\right) = -\frac{d}{dx}\left(-6x + 18\right) = -(-6) = 6$$

3. **Combine Radial Terms:**
   $$R_2^1(r) = \frac{1}{12\sqrt{6}\,a_0^{3/2}} \times \left(\frac{r}{a_0}e^{-r/(2a_0)}\right) \times 6 = \frac{1}{2\sqrt{6}}\,a_0^{-3/2}\left(\frac{r}{a_0}\right)e^{-r/(2a_0)}$$

#### Step 2: Angular Part $Y_1^1(\theta,\phi)$
$$Y_1^1(\theta,\phi) = (-1)^1 \sqrt{\frac{(2(1)+1)(1-1)!}{4\pi(1+1)!}}\,P_1^1(\cos\theta)e^{i\phi} = -\sqrt{\frac{3(0!)}{4\pi(2!)}}\,P_1^1(\cos\theta)e^{i\phi} = -\frac{1}{2}\sqrt{\frac{3}{2\pi}}\,P_1^1(\cos\theta)e^{i\phi}$$

1. **Associated Legendre Polynomial $P_1^1(\cos\theta)$:**
   $$P_1(x) = \frac{(-1)^1}{2^1(1!)}\frac{d}{dx}(1 - x^2) = -\frac{1}{2}(-2x) = x$$
   $$P_1^1(x) = (-1)^1 (1 - x^2)^{1/2}\frac{d}{dx}(x) = -\sqrt{1 - x^2}$$
   Substituting $x = \cos\theta$:
   $$P_1^1(\cos\theta) = -\sqrt{1 - \cos^2\theta} = -\sin\theta$$

2. **Combine Angular Terms:**
   $$Y_1^1(\theta,\phi) = \left(-\frac{1}{2}\sqrt{\frac{3}{2\pi}}\right)(-\sin\theta)e^{i\phi} = +\frac{1}{2}\sqrt{\frac{3}{2\pi}}\sin\theta\,e^{i\phi}$$

#### Step 3: Final Wave Function $\psi_{211}(r,\theta,\phi)$
Multiplying $R_2^1(r)$ and $Y_1^1(\theta,\phi)$ together:
$$\psi_{211}(r,\theta,\phi) = \left[\frac{1}{2\sqrt{6}}\,a_0^{-3/2}e^{-r/(2a_0)}\left(\frac{r}{a_0}\right)\right] \left[\frac{1}{2}\sqrt{\frac{3}{2\pi}}\sin\theta\,e^{i\phi}\right]$$
$$\psi_{211}(r,\theta,\phi) = \frac{1}{4}\sqrt{\frac{3}{12\pi}}\,a_0^{-3/2}e^{-r/(2a_0)}\left(\frac{r}{a_0}\right)\sin\theta\,e^{i\phi} = \frac{1}{8\sqrt{\pi}}\,a_0^{-3/2}e^{-r/(2a_0)}\left(\frac{r}{a_0}\right)\sin\theta\,e^{i\phi}$$