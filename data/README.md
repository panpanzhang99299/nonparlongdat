# Mock Data Description

The mock data is generated from a linear mixed effects model

\[Y_{ij} = \theta_0 + \theta_1 X_{i} + \theta_2 Z_{ij} + b_i + \varepsilon_{ij}.\]

The parameters are specified with
- $\theta_1 = 12$
- $\theta_2 = 3$
- $\theta_3 = -0.3$
- $\sigma_b = 0.8$
- $\sigma_{\varepsilon} = 1.2$ 

There are a total of 450 subjects in the mock data, and each has 3 visits.

The missing data is generated under an MAR mechanism.