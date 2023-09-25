# radialvelocity
The radial velocity technique is the cornerstone of modern astronomy. Before the advent of NASA's *Kepler* efforts, radial velocity measurements served as the primary method of exoplanet detection: **Campbell et al. 1988, Latham et al. 1989, Mayor & Queloz 1995, Marcy & Butler 1996., Akeson et al. 2013.** The mathematical solution of the two-body problem shows that when a planet orbits a star in an elliptical path, the star also moves in an ellipse about a common center of mass. The technique relies on measuring the Doppler shifts in the line spectra of stars and assessing their recessional velocities in these elliptical orbits. This is calculated using the relativistic Doppler equation.
```math
\begin{equation}
    z \equiv \frac{\lambda_{\mathrm{obs}} - \lambda_{\mathrm{rest}}}{\lambda_{\mathrm{rest}}} = \frac{1}{\gamma(1+v/c)} - 1 \quad \mathrm{where} \quad \gamma = \frac{1}{1-(v/c)^2} - 1 \tag{1.1}
\end{equation}
```
where $\gamma$ is the Lorentz factor. Then fitting this data for the orbital parameters yields information of any bodies orbiting it. The velocity, $v$, is a function of time which takes a periodic solution since for one half of the orbit, the spectral lines will be blueshifted and for the other half, redshifted. The motion of a body around another can then be understood in the context of the two-body problem; the solution to which was first presented by Sir Isaac Newton in 1687 [1]. Newton demonstrated that the elliptical motion of bodies and Kepler's three laws (1609, 1619) could both be described by his law for gravitation. But what Newton achieved was only showing that Kepler's observations could be explained by an inverse-square law. To actually find the positions and velocities of bodies in the two-body problem is not quite as straightforward; this is termed the *Kepler* problem. 

The Keplerian problem is six-dimensional, describing orbits in a binary system by six parameters. These are shown in the table below.

| Quantity                    | Parameter  | Definition                                                                                   |
|-----------------------------|------------|----------------------------------------------------------------------------------------------|
| Eccentricity                | $e$        | The ratio of the distances from the center of the ellipse to one of the foci and to one of the vertices of the ellipse              |
| Period                      | $P$        | The time taken for a body to complete one revolution                                           |
| Inclination                 | $i$        | The angle between an orbital plane and a reference plane                                      |
| Ascending node              | $\Omega$   | The angle made between the reference direction and the point where the orbiting body crosses the reference plane going "upward" |
| Argument of Periastron      | $\omega$   | The angle between the ascending node and the periapsis.                                       |
| Epoch of periastron         | $T_0$      | A moment in time used as a reference point for time-varying astronomical quantities.|

These parameters can completely define an orbit but are not yet whole. The *Keplerian* problem can be reduced to the solutions of primarily three equations, [2].

```math
\begin{gather}
    M = E - e\sin E \tag{1.2}\\
    \nu = 2\arctan \left(\sqrt{\frac{1+e}{1-e}}\tan\frac{E}{2}\right) \tag{1.3} \\
    \dot{\mathrm{z}} = \mathrm{v} = \kappa + K\left[\cos(\nu(t) + \omega) + e\cos(\omega)\right] \tag{1.4}
\end{gather}
```
where $M$ is the mean anomaly, $E$ is the eccentric anomaly, $\nu$ is the true anomaly, $K$ is the velocity semi-amplitude, i.e., a measure of the "wobble" induced on the star from an orbiting body,
```math
K = \frac{n a\sin i}{\sqrt{1-e^2}},
```
$\kappa$ is the long-term mean velocity from the binary, $a$ is the semi-major axis, and $n$ is the mean motion $n = 2\pi/P$. Equation 1.2 is known as *Kepler's* equation and no closed-form solution to it exists. Physically, Kepler's equation demonstrates the relationship between the polar coordinates defining a celestial body and the time elapsed from a given moment. Equation 1.3 is derived from arguments of the true anomaly which goes beyond the mathematical scope of this paper [5]. Equation 1.4 is the radial velocity equation. The data is fit for the parameters in that equation. All the parameters outlined in Table 1 are shown in a diagram of the Keplerian orbit below. It is worth noting that $\Omega$ does not appear in expressions that solve for the radial velocity. This is obvious since the orientation of the orbit in the plane of the sky does not affect the observations. This also extends itself to imply that $\Omega$ cannot be found from radial velocity data.

![image](https://github.com/yashpincha/radialvelocity/assets/142161252/f81cfd95-19e2-4e39-aec2-d504b1e72f69)

**Kepler Orbital Parameters**
*Left:* A three-dimensional view of the Kepler problem with relevant parameters labeled. The apastron is the point of furthest distance from the center of gravity. The reference plane is considered to be equatorial, i.e., when the inclination is equal to zero. This is often a very difficult quantity to measure, and therefore only the projected value of the semi-major axis $a\sin i$ is measurable. $\gamma$ denotes the direction pointing towards the observer, i.e., Earth. The diagram is not drawn to scale. *Right:* A two-dimensional projection of the Kepler problem with the **true anomaly** indicated. The inferior conjunction is the closest point at which the planet crosses the observer's line of sight.

As the star moves along the smaller elliptical orbit depicted by the dashed line in the Left figure, its atomic spectra will shift as well. The radial velocity of celestial bodies are measured by taking high-resolution spectra; this is typically done by superimposing Echelle spectograph. These are then compared with laboratory-measured spectral line wavelengths. Radial velocity measurements made in this manner are termed as **absolute** radial velocities. Such measurements are limited by the wavelength calibration of the spectrograph and typically have accuracies of the order of 100m/s (Chubak et al., 2012). An alternative method is to measure differential radial velocities, i.e., the change in redshift between two epochs. Another factor worth considering is that as observatory taking the spectographs is itself moving --- the Earth rotates at ~300m/s and orbits the Sun at ~30km/s --- the measured radial velocity of a perfectly stable star will appear to vary on diurnal and annual timescales. Since even the Earth takes an elliptic orbit, its velocity is not constant and would increase at the point of perihelion. Such measurements of radial velocity are corrected such that the measurements are as an observer would have measured in an inertial reference frame. This is termed as *barycentric* correction and the method of correcting it follows from Wright and Eastman (2014).

## Orbital anomalies 
<p align="center">
  <img width="440" height="380" src="https://github.com/yashpincha/radialvelocity/assets/142161252/a6063adf-2605-48cd-a958-d204cdde19a7">
</p>
The anomalies are angular measures that define the position of a mass moving _along_ a Keplerian orbit. These are derived from the seven primary parameters in Table 1. The mean anomaly, $M$, is the fraction of an elliptical orbit's period that has elapsed since the orbiting body passed periastron. It is the angular distance from the center of gravity which a fictitious body would have if it moved in a circular orbit, with constant speed, in the same orbital period as the actual body in its elliptical orbit [3, 4]. Kepler’s second law postulates that planets sweep equal areas in equal times --- that is to say that celestial bodies do not move uniformly in orbit but speed up at periastron and slow down at the apastron. To circumvent this irregularity, the mean anomaly is introduced. It is far easier to calculate angles in a circular orbit and later using it with the eccentric anomaly to arrive at the true anomaly. The eccentric anomaly, $E$, is one of the angles on the right-angled triangle where one vertex lies at the center of the ellipse, the adjacent side lies on the semi-major axis, the hypotenuse of length $a$ lies on the auxiliary circle of the ellipse, and the opposite side is perpendicular to the semi-major axis passing through the point $P$, the location of the orbiting body. It is essentially like the true anomaly but with the exception of the center being located at the center of the ellipse as opposed to the center of mass. The true anomaly is the angle made between the semi-major axis and the position of the orbiting body as seen from the periastron.

## Mass determination
The semi-amplitude, $K$, and the period, $P$, are both required in determining the mass of the exoplanet. Writing out Kepler's third law,
```math
    G(M_1+M_2) = \omega^2 a^3
```
where $M_1$ is the mass of the star and $M_2$ is that of the unseen body, i.e., the exoplanet, $a$ is equal to the two semi-major axes, $a = a_1 + a_2$, G is the gravitational constant, G = $6.67\times 10^{-11} \mathrm{m}^3\mathrm{kg}^{-1}\mathrm{s}^{-2}$, $\omega$ is the orbital frequency. Using the fact that the binary system's center of mass satisfies
```math
    M_1a_1 = M_2a_2
```
$a$ can be rewritten as 
```math   a = a_1 + a_2 = a_1\qty(1 + \frac{a_2}{a_1}) = a_1\qty(1 + \frac{M_1}{M_2}) = \frac{a_1(M_1+M_2)}{M_2}
```
which when substituted in Kepler's third law yields
```math
    G(M_1+M_2) = \omega^2\left(\frac{a_1(M_1+M_2)}{M_2}\right)^3
```
which is the same as
```math
    {\displaystyle {\frac {M_{2}^{3}}{M_{\mathrm {tot} }^{2}}}={\frac {\omega^{2}a_{1}^{3}}{G}}.}
```
Since the Doppler method only measures the $\sin$ component of the recessional/processional velocities, the measurements depend on the orbital inclination. When $i = 0$, the orbit is being viewed face-on and the sine component of the velocity is zero. When $i=90^\circ$, the orbit is being viewed edge-on. The peak velocity is the same as the semi-amplitude of the curve. For a circular orbit, when eccentricity is zero, the peak velocity is given by
```math
    K_1 = \mathrm{v_1}\sin i = \omega a_1\sin i
```
Substituting this for $a_1$ yields
```math
    {\displaystyle {\frac {M_{2}^{3}}{M_{\mathrm {tot} }^{2}}}={\frac {K^{3}}{G\omega\mathrm {sin} ^{3}i}}.}
```
The mass function, then, defined as $f(M_1, M_2)$ is
```math
    f=\frac {(M_2\sin i)^3}{(M_{1}+M_{2})^{2}}=\frac {PK^{3}}{2\pi G} = \frac{4\pi^2 (a_1\sin i)^3}{GP^2}
```
For a known mass $M_1$, an estimate for $M_{2, {\mathrm{min}}}$ can be determined. The true mass $M_2$ depends on the orbital inclination. This is an inherent disadvantage of spectroscopic methods of studying binaries but a lower bound for the mass still offers sufficient data for study. The inclination can typically be determined from eclipses or be modelled using ellipsoidal variations [7]. **Limits.** In the case of $M_1 \gg M_2$, as is the case with exoplanets, the mass function approximates to 
```math
    {\displaystyle f\approx {\frac {M_{2}^{3}\ \mathrm {sin} ^{3}i}{M_{1}^{2}}}.}
```
And in the case of $e \neq 0$, the mass function is [8],
```math
    f \equiv \frac{(M_2\sin i)^3}{(M_1 + M_2)^2} = \frac{PK^3}{2\pi G}(1-e^2)^{(3/2)}
```
the equations outlined in this section are used to determine the masses of orbiting exoplanets.

## Lomb-Scargle
To fit the raw data, an initial guess of the period is required. To achieve this, the Lomb-Scargle periodogram is used. The Lomb-Scargle periodogram is a popular method used for detecting periodic signals in unevenly-sampled data. Alternative methods which are generally more accurate include Bayesian approaches, Fourier methods, and phase-folding methods. The LS periodogram is far more convenient and is a least-squares method - a least-squares method is a process of finding the best-fitting curve by reducing the sum of the squares of the offsets. The LS periodogram primarily works on this principle. Its mathematical formulation is as directed in [9].

## How we actually get there
To determine the orbital parameters, the raw radial velocity data is fit to a cosine of the earlier form, minimizing the $\chi^2$. The initial guess for $\kappa$ is first made by averaging the radial velocity data: $\kappa = \sum V_i/N$ where N is the number of data points. Some straightforward mathematical manipulation yields the other parameters. From the cosine4, it is obvious that the radial velocity is a maximum when $\cos(\nu + \omega) = 1$ and a minimum when $\cos(\nu + \omega) = -1$. At the maximum and minimum radial velocities, the equation reduces to 

```math
\frac{V_{\mathrm{max}} + V_{\mathrm{min}}}{2} = \kappa + Ke\cos\omega
```
which ultimately yields
```math
        E_\mathrm{max} = 2\arctan\left[\left(\frac{1+e}{1-e}\right)^{-1/2}\tan\left(\frac{-\omega}{2}\right)\right]
```
## Exoplanet candidate criteria

A literature measurement for period $P$, eccentricity $e$, the semi-amplitude $K$, and stellar mass $M_\star$ should exist. A value for the planetary mass should exist, $m\sin i$. Sufficient ancillary radial velocity data, $N > 40$, should exist. The data should have been updated no later than 2000. The exoplanet must orbit a K-type star (not particularly relavant; this was just the paper's focus since K-type stars are most likely to hold habitable planets).

| Identity       | Source                | $M_\star$ ($M_\odot$) | $\Delta M_\star$ ($M_\odot$) |  ${M_{\star,\mathrm{avg}}}$ ($M_\odot$) | $\Delta_t M_\star$ ($M_\odot$) | Velocity           |
|----------------|-----------------------|------------------------|------------------------------|-----------------------------------------|-------------------------------|---------------------|
| 6 Lyncis       | Luhn et al. [16]      | 1.44                   | 0.00                         | 1.61                                    | 0.18                          | [35]                |
|                | Sousa et al. [17]      | 1.56                   | 0.05                         |                                         |                               |                     |
|                | Bowler et al. [18]     | 1.82                   | 0.13                         |                                         |                               |                     |
| 11 Ursae Minoris | Baines et al. [21]  | 1.43                   | 0.22                         | 1.92                                    | 0.45                          | [19]                |
|                | Dollinger et al. [19] | 1.80                   | 0.25                         |                                         |                               |                     |
| 14 Andromedae  | Sato et al. [20]      | 1.12                   | 0.24                         | 1.16                                    | 0.43                          | [35]                |
|                | Sousa et al. [17]      | 1.20                   | 0.19                         |                                         |                               |                     |
| BD+14 4559       | TICv8 [22]            | 0.79                   | 0.09                         | 0.82                                    | 0.35                          | [24]                |
|                  | Niedzielski et al. [24]| 0.82                   | 0.11                         |                                         |                               |                     |
|                  | Santos et al. [23]     | 0.86                   | 0.15                         |                                         |                               |                     |
| BD-10 3166       | TICv8 [22]            | 0.92                   | 0.11                         | 0.92                                    | 0.16                          | [36]                |
|                  | Takeda et al. [25]     | 0.92                   | 0.05                         |                                         |                               |                     |
| Epsilon Borealis | Mortier et al. [26]  | 1.44                   | 0.18                         | 1.57                                    | 0.28                          | [27]                |
|                  | Lee et al. [27]        | 1.70                   | 0.11                         |                                         |                               |                     |
| Epsilon Eridani | Llop-Sayson et al. [28] | 0.82                 | 0.02                         | 0.83                                    | 0.10                          | [36]                |
|                  | Rosenthal et al. [29]  | 0.82                   | 0.02                         |                                         |                               |                     |
|                  | Takeda et al. [25]     | 0.85                   | 0.00                         |                                         |                               |                     |
|                  | Benedict et al. [30]   | 0.83                   | 0.05                         |                                         |                               |                     |
| HD 2952          | Stassun et al. [39]    | 1.97                   | 0.38                         | ~                                       | ~                             | [39]                |
| HD 3651           | Wittenmyer et al. [31] | 0.88                   | 0.03                         | 0.88                                    | 0.10                          | [36]                |
|                  | Takeda et al. [25]     | 0.88                   | 0.03                         |                                         |                               |                     |
|                  | Rosenthal et al.       | 0.88                   | 0.04                         |                                         |                               |                     |
| HD 32518         | Sousa et al. [17]      | 1.16                   | 0.16                         | 1.15                                    | 0.34                          | [19]                |
|                  | Dollinger et al. [19] | 1.13                   | 0.18                         |                                         |                               |                     |
| HD 96127        | Sousa et al. [17]      | 1.29                   | 0.27                         | 1.10                                    | 0.52                          | [32]                |
|                  | Gettel et al. [32]     | 0.91                   | 0.25                         |                                         |                               |                     |
| HD 130322        | Hinkel et al.          | 0.92                   | 0.03                         | 0.92                                    | 0.27                          | [37]                |
|                   | TICv8 [22]            | 0.93                   | 0.13                         |                                         |                               |                     |
|                  | Ghezzi et al.          | 0.90                   | 0.11                         |                                         |                               |                     |
| HD 192263         | Rosenthal et al. [29]  | 0.82                   | 0.03                         | 0.81                                    | 0.25                          | [23]                |
|                 | TICv8 [22]            | 0.82                   | 0.10                         |                                         |                               |                     |
|                 | Dragomir et al. [33]   | 0.81                   | 0.02                         |                                         |                               |                     |
| HD 240237       | Sousa et al.           | 0.61                   | 0.08                         | 1.15                                    | 0.49                          | [32]                |
|                 | Gettel et al. [32]     | 1.69                   | 0.41                         |                                         |                               |                     |

## Results for some exoplanets. 

| Identity       | K       | ΔK    | Period   | ΔP    | e     | Δe   | a     | Δa   | $m\sin i$ (Jupiter masses) | Δm    | $\chi^2$   |
|----------------|---------|-------|----------|-------|-------|------|-------|------|-----------------------|-------|------------|
| 6lyn b         | 36.206  | 1.532 | 899.229  | 13.971 | 0.135 | 0.044| 2.116 | 0.218| 2.296                 | 0.480 | 62.124     |
| 11UMi b        | 192.258 | 1.117 | 516.922  | 0.662  | 0.081 | 0.006| 1.480 | 0.143| 10.425                | 2.021 | 1534.049   |
| 14And b        | 100.028 | 1.465 | 185.842  | 0.239  | 0.008 | 0.015| 0.669 | 0.077| 3.097                 | 0.715 | 331.191    |
| BD+14 4559     | 50.329  | 1.694 | 264.572  | 0.880  | 0.329 | 0.028| 0.756 | 0.108| 1.319                 | 0.380 | 164.050    |
| BD-10 3166     | 59.770  | 2.074 | 3.488    | 0.000  | 0.018 | 0.016| 0.044 | 0.004| 0.430                 | 0.079 | 76.435     |
| Epsilon CrB    | 134.609 | 2.530 | 420.648  | 0.857  | 0.069 | 0.017| 1.277 | 0.076| 6.688                 | 0.805 | 197.310    |
| HD 2952 b      | 272.277 | 0.841 | 272.277  | 0.421  | 0.229 | 0.034| 0.066 | 0.066| 0.971                 | 0.132 | 1262.922   |
| HD 3651 b      | 16.210  | 0.432 | 62.192   | 0.005  | 0.672 | 0.013| 0.295 | 0.010| 0.215                 | 0.215 | 2588.844   |
| HD 32518 b     | 110.330 | 1.985 | 157.567  | 0.202  | 0.048 | 0.017| 0.597 | 0.059| 3.208                 | 0.637 | 210.766    |
| HD 96127 b     | 92.784  | 1.668 | 582.390  | 1.217  | 0.560 | 0.012| 1.409 | 0.222| 3.366                 | 1.063 | 3683.564   |
| HD 130322 b    | 115.634 | 1.535 | 10.724   | 0.005  | 0.049 | 0.013| 0.092 | 0.009| 1.183                 | 0.224 | 185.637    |
| HD 192263 b    | 60.215  | 0.698 | 24.347   | 0.003  | 0.018 | 0.012| 0.153 | 0.015| 0.747                 | 0.151 | 653.077    |
| HD 240237 b    | 90.121  | 3.270 | 746.960  | 2.821  | 0.396 | 0.028| 1.689 | 0.240| 4.061                 | 1.162 | 1142.857   |

## Directions
Find some csv files that I could salvage in the `csv` folder. The citations are in `cite.txt`. The `images` folder contains an example of expected plots. The results aren't the *most* accurate yet but the MCMC algorithm does a much better job than the least sq. one, as expected.

Yash Pincha, 2022. All images by author.
