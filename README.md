# NCrystal plugin Liquids

This is a preliminary plugin for the inclusion of liquid modelling on-the-fly in NCrystal. Example data files are included in the data directory.

To use the plugin, a `@CUSTOM_LIQUIDS` section needs to be appended to the input of a `.ncmat` file. An example of the parameter options are given below.
```
@CUSTOM_LIQUIDS
<atom>  <model>  <c>  <wt>  <ws>
yk_model  <state>          # <state> = ortho | para  (liquid H2/D2 only)
<Q_1>   <S_1>              # S(Q) rows; omit when <model> = incoherent
  ...
<Q_n>   <S_n>
```

*The `@CUSTOM_LIQUIDS` section. Angle-bracketed terms are user inputs.*

Each entry begins with a single line specifying the atom and its scattering parameters, optionally followed by a `yk_model` line, which selects the Young–Koppel model for liquid hydrogen and deuterium, and a table of structure-factor values. The block is repeated for each atom (e.g. `D` then `O` for D₂O). The variables are:

- `<atom>`: atom the entry applies to (`H`, `D`, `O`, …);
- `<model>`: incoherent/coherent model selection, can be either `incoherent`, `skold`, or `vineyard`;
- `<c>`: diffusion constant;
- `<wt>`: translational weight;
- `<ws>`: solid-like weight;
- `<state>`: Young–Koppel rotational state, either `ortho` or `para`; used only for liquid H₂/D₂;
- `<Q_i>`, `<S_i>`: structure-factor rows, with *Q* in Å⁻¹ and *S(Q)* dimensionless; omitted when `<model>` is `incoherent`.

# Reference
D. D. DiJulio, J. I. Marquez Damian, T. Kitelmann and G. Muhrer, A liquids plugin for NCrystal. Presented at ICANS XXV  (13-17 April 2026). https://indico.ess.eu/event/3811/contributions/24554/
