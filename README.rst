# MF-DXA Multifractal Detrended Cross-Correlation Analysis

C++ implementation of Multifractal detrended cross-correlation analysis for two nonstationary signals from [1].

[1] Zhou, W. X. (2008). "Multifractal Detrended Cross-Correlation Analysis for Two Nonstationary Signals". *Physical Review E—Statistical, Nonlinear, and Soft Matter Physics*, 77(6), 066211.

## 1. Data Preparation

- Start with two time series $x(i)$ and $y(i)$, where $i = 1, 2, \ldots, M$.
- Ensure that both series have the same length $M$.

## 2. Calculation of Profiles

- Calculate the profiles $X(i)$ and $Y(i)$ by subtracting the means from the original series and cumulatively summing:
  $$
  X(i) = \sum_{k=1}^i [x(k) - \bar{x}]
  $$
  $$
  Y(i) = \sum_{k=1}^i [y(k) - \bar{y}]
  $$
- Where $\bar{x}$ and $\bar{y}$ are the means of $x$ and $y$, respectively.

## 3. Segmentation

- Divide both series into $N_s = \text{int}(N/s)$ non-overlapping segments of length $s$.

## 4. Calculation of Local Covariance

- For each segment $v$, calculate the local covariance $F^2_{xy}(s,v)$ using polynomials of order $m$ (typically linear or quadratic) to detrend both series.
- The general formula is:
  $$F^2_{xy}(s,v) = \frac{1}{s} \sum_{k=1}^s \left[ \left(X((v-1)s+k) - \tilde{X}_v(k)\right) \cdot \left(Y((v-1)s+k) - \tilde{Y}_v(k)\right) \right]$$
- Where $\tilde{X}_v(k)$ and $\tilde{Y}_v(k)$ are the estimated local trends.

## 5. Calculation of Fluctuation Function

- Calculate the $q$-th order fluctuation function:
  $$
  F_q(s) = \left( \frac{1}{N_s} \sum_{v=1}^{N_s} \text{sign}\left[F_{xy}^2(s,v)\right] \cdot \left| F_{xy}^2(s,v) \right|^{\frac{q}{2}} \right)^{\frac{1}{q}}
  $$  
- Repeat this calculation for different values of $q$ (both positive and negative).

## 6. Repetition for Different Scales

- Repeat steps 3-5 for different segment lengths $s$.

## 7. Scaling Law Analysis

- Plot $\log(F_q(s))$ vs $\log(s)$ for each value of $q$.
- If the series are cross-correlated in a multiscale manner, you should observe a linear relationship.
- The slope of this linear relationship is the generalized Hurst exponent $h(q)$.

## 8. Multiscale Spectrum Calculation

- Calculate the multiscale spectrum $\tau(q) = qh(q) - 1$.
- Calculate the generalized fractal dimension $D(q) = \frac{\tau(q)}{q - 1}$.

## 9. Result Interpretation

- Analyze how $h(q)$, $\tau(q)$, and $D(q)$ vary with $q$ to characterize the multiscale properties of the cross-correlation between the two series.

## Notes

- For very long time series, computational optimization techniques may be necessary.
- The choice of polynomial order for detrending can influence the results, so experimenting with different orders may be required.
- It is important to check the robustness of the results by varying the analysis parameters, such as the range of $s$ and $q$.

# HOW TO USE

### Steps to Install Armadillo and Set Up a Project

1. **Install Armadillo:**
   - Download Armadillo from the official website: [https://arma.sourceforge.net/](https://arma.sourceforge.net/).
   - Follow the installation instructions provided on the Armadillo website to properly install the library on your system.

2. **Modify the Makefile to Include Armadillo:**
   - Open the `Makefile` for your project.
   - Update the library and include paths to point to the locations where Armadillo is installed.
     - Example:
       ```make
       IMPLIB = -L/path/to/armadillo/lib
       IMPINC = -I./include -I/path/to/armadillo/include
       ```

3. **Change the `examples/test.cpp` Path in the Makefile:**
   - Locate the part of the `Makefile` that references `examples/test.cpp`.
   - Replace this path with the path to your own script (e.g., `your_script.cpp`):
     ```make
     SOURCE = path/to/script.cpp src/mfdxa.cpp
     ```

4. **Compile Using Make:**
   - Once the `Makefile` is properly set up with the correct paths, compile your program by running:
     ```bash
     make
     ```

5. **Plot a Graph Using Gnuplot:**
   - To visualize your results, you can use Gnuplot. After running your program, you can pipe or redirect the output to Gnuplot for plotting (see /examples)
