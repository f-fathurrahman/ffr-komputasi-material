## The second slide

These bullets will be displayed incrementally.

::: {.incremental}

- Get in bed
- Count sheep

:::

Equation: $\alpha + \beta$

## Schroedinger equation

::: {.fragment .fade-in}
$$
\mathrm{i}\hbar \frac{\partial}{\partial t} \left| \psi(x,t) \right\rangle =
\hat{H} \left| \psi(x,t) \right\rangle
$$
:::

::: {.fragment .fade-in}
Colored equation, using non-standard LaTeX, MathJaX version 2.
$$
\hat{H} = 
\color{green}{-\frac{\hbar^2}{2m} \nabla^2} +
\color{red}{V(x,t)}
$$ 
:::



## Example underbrace {.smaller}

$$
\hat{H} = 
\underbrace{\color{green}{-\frac{\hbar^2}{2m} \nabla^2}}_{\text{kinetic}} +
\underbrace{\color{red}{V(x,t)}}_{\text{potential}}
$$ 


## Example cases

$$
 u(x) = 
  \begin{cases} 
   \exp{x} & \text{if } x \geq 0 \\
   1       & \text{if } x < 0
  \end{cases}
$$


## Example arrow

$$
 A \xleftarrow{\text{this way}} B 
  \xrightarrow[\text{or that way}]{ } C
$$



## Using HTML?

<p style="color: red;">Some paragraph</p>

$$
{\color{blue}\alpha} + \beta
$$


## Fragment

::: {.fragment fragment-index=3}
Appears last
:::

::: {.fragment fragment-index=1}
Appears first
:::

::: {.fragment fragment-index=2}
Appears second
:::



## A code

Example Julia code:

```julia
function hello()
    println("Hello from efefer")
end
```




## LaTeX Equations

[MathJax](https://www.mathjax.org/) rendering of equations to HTML

::: columns
::: {.column width="40%"}
``` tex
\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}
\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}
```
:::

::: {.column width="60%"}
```{=tex}
\begin{gather*}
a_1=b_1+c_1\\
a_2=b_2+c_2-d_2+e_2
\end{gather*}
```
```{=tex}
\begin{align}
a_{11}& =b_{11}&
  a_{12}& =b_{12}\\
a_{21}& =b_{21}&
  a_{22}& =b_{22}+c_{22}
\end{align}
```
:::
:::

::: footer
Learn more: [LaTeX Equations](https://quarto.org/docs/authoring/markdown-basics.html#equations)
:::

## Column Layout {.smaller}

Arrange content into columns of varying widths:

::: columns
::: {.column width="35%"}
#### Motor Trend Car Road Tests

The data was extracted from the 1974 Motor Trend US magazine, and
comprises fuel consumption and 10 aspects of automobile design
and performance for 32 automobiles.
:::

::: {.column width="3%"}
:::

::: {.column width="62%"}
```r
knitr::kable(head(mtcars)[,c("mpg",	"cyl", "disp", "hp", "wt")])
```
:::
:::

::: footer
Learn more: [Multiple Columns](https://quarto.org/docs/presentations/revealjs/#multiple-columns)
:::


## Line Highlighting {.smaller}

-   Highlight specific lines for emphasis
-   Incrementally highlight additional lines

``` {.python code-line-numbers="4-5|7|10"}
import numpy as np
import matplotlib.pyplot as plt
r = np.arange(0, 2, 0.01)
theta = 2 * np.pi * r
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r)
ax.set_rticks([0.5, 1, 1.5, 2])
ax.grid(True)
plt.show()
```

:::{.incremental}
- Explanation 1
- Explanation 2
- Explanation 3
:::

::: footer
Learn more: [Line Highlighting](https://quarto.org/docs/presentations/revealjs/#line-highlighting)
:::