# Dammage functions:
The authors used differents definition for calculating environment damages. This section presents the diferent functions the user can use to calculate them, depending of the method he prefers.

## ACE
The damage function in ACE takes an exponential functional form: 
``D(T_{1,t})=1-\exp(-\xi_0\exp(\xi_1 T_{1,t})+\xi_0)``
with ``\xi_0``, a free damage parameters and the paramteter ``\xi_1=\frac{\log 2}{s}`` pinned down by climate sensitivity *s*.

Input:
- ``\xi_0``: xi0 (Calibrated to 0.022 in the paper. The base calibration is an exact match of the two calibration points 0° and 2.5° in the 2007 model)
- ``\xi_1``: xi1 (Estimated to ``\frac{1}{4}`` in the paper)
- ``T``: temperature

Syntax:
```julia
dam_ACE(xi0, xi1, T)
```
Output:
The calculated ACE damage estimate.

## DICE
DICE damage function replace the quadratic damage term ``aT^2``, to limit damages to 100 percent of production. 
``D(T)=1-\frac{1}{1+aT^2}``

Input:
- ``a`` : Damage coefficients (``a=0.0028`` for 2007, ``a=0.00267`` for 2013, ``a=0.00236`` for 2016)
- ``T``: temperature

Syntax:
```julia
ACE_Traeger_replication.dam_DICE(a, T)
```
Output:
The calculated DICE damage estimate.

## Howard and Sterner (2017)
This damage function is based on Howard and Sterner (2017):
``D(T)=aT^2``

Input:
- ``a`` : Damage coefficient. By default, it is set to 0.01145. Write it explicitely, when you are using the dam_Sterner function if you want another ``a``.
- ``T``: temperature

Syntax:
```julia
ACE_Traeger_replication.dam_Sterner(T, a)
```

Output:
The calculated DICE-Howard-Sterner damage estimate.


The authors used differents 