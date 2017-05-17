# Line-VISAR analysis routine

## Theory

The fringe shift recorded is directly proportional to the velocity of the target. It is, therefore, the goal of the analysis routine to extract the percentage fringe shift (at a given location and time). Several methods exist for accomplishing this. However, the Fourier transform method (described below), first proposed by Takeda et al. \cite{M.-Takeda:1983bh}, has been determined to be the most accurate \cite{Celliers:2004uo}. 

The raw data (fig. \ref{fig:raw}) is typically cropped to analyze the region of planar shock breakout (Figure \ref{fig:visar_analy}-a). This image can be described mathematically using equation \ref{eqn:2dgen}, with $ b(x,t ) = I_1+I_2 $ being the background and $ a(x,t) =  2\mathbf{E_1}\mathbf{E_2} $ describing the intensity of the fringes. $\phi (x,t)$ represents the phase of the fringes and $2\pi f_{0}x + \delta_{0} $ describes the linear phase ramp of the background fringe pattern. The goal is to find $ \phi (x,t) $, which is directly proportional to the velocity of the target via equation \ref{eqn:fshift} and equation \ref{eqn:2dvpf}. Rewriting equation \ref{eqn:2dgen} in terms of its complex components yields:

$$
\begin{align}
&& f(x,t) &= b(x,t)+c(x,t)e^{i2\pi f_{0}x} + c^*(x,t)e^{-i2\pi f_{0}x },& \\
\text{with} && c(x,t) &= \frac{1}{2}a(x,t)e^{i\delta_{0}}e^{i\phi (x,t)}. 
\end{align}
$$

Applying a Fourier transform to the data at each point-in-time (Figure \ref{fig:visar_analy}-b) allows for the filtering of specific frequencies, such that the background, $b(x,t)$, can be removed by setting the pixel values to zero, as illustrated in Figure \ref{fig:visar_analy}-c and  \ref{fig:visar_analy}-d:

$$
\begin{align*}
F(f,t) &= B(f,t)+\int_{ -\infty }^{\infty} c(x,t) e^{i2\pi f_{0}x}e^{-ifx} \; dx + \int_{-\infty }^{\infty}c^*(x,t)e^{-i2\pi f_{0}x }e^{-ifx} \; dx \\[1.5ex]
&=B(f,t)+\int_{-\infty }^{\infty}c(x,t)e^{i2\pi (f_{0}-f) x} \; dx + \int_{-\infty }^{\infty}c^*(x,t)e^{-i2\pi (f_{0} +f )x} \;dx \\[1.5ex]
&=\cancelto{0}{B(f,t)}+\cancelto{0}{C^*(f+f_0,t)}+C(f-f_0,t).
\end{align*}
$$

Applying an inverse Fourier transform (Figure  \ref{fig:visar_analy}-e), this can be written as:

$$
\begin{align}
d(x,t) &= \int_{-\infty}^{\infty}C(f-f_0,t)e^{ixf} \; df \nonumber \\
&= \int_{-\infty}^{\infty}C(f-f_0,t) \left ( \cos(xf) + i\sin(xf) \right ) \; df  \label{eqn:ifft} \\
&=c(x,t)e^{2\pi i f_0x} \label{eqn:filtered}
\end{align}
$$

Equation \ref{eqn:filtered} has both a real and imaginary valued functions (as seen in eqn. \ref{eqn:ifft}). Where

$$
\begin{align}
&& \operatorname{Re} [d(x,t) ] &\propto  \sin( \phi (x,t) + 2\pi f_0 x + \delta_0) \label{eqn:re},& \\ 
\text{and} && \operatorname{Im} [d(x,t) ] &\propto  \cos( \phi (x,t) + 2\pi f_0x + \delta_0) \label{eqn:im},
\end{align}
$$

are $\pi/2$ out of phase. Taking the $\arctan$ of the ratio (Figure  \ref{fig:visar_analy}-f) allows the phase, $\phi (x,t) + 2\pi f_0x + \delta_0$, to be extracted: 

$$
\begin{equation*}
W( \phi (x,t) + 2\pi f_0x + \delta_0) = \arctan\left ( \frac{\operatorname{Re} [d(x,t) ]}{\operatorname{Im} [d(x,t) ]} \right ).
\end{equation*}
$$

The resulting function $W$ has discontinuities representing $\pi$ shifts as the $\arctan$ moves through full rotations. The velocity signal can be constructed by setting these shifts to zero and scaling the values by the proportionality factor VPF (eqn: \ref{eqn:2dvpf}). The programatic method for reconstructing the velocity trace from the time dependent values in wrapped\_phase parameters can be accomplished via the (Python) script bellow:
\vspace{4mm}

## Implementation of Analysis Routine

