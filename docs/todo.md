# To do

- Handle the 4-angle case completely, both in simulation and reduction. (Currently, it is implemented in the reduction, albeit it not tested, since no actual or simulated data are available.)

- Implement a 2-angle (polarizer) case.

- Include the configuration files with the package, as constant strings in a module that can be printed by the user.

- Simplify interpolation correction (remove AstroAlign dependency; restrict to shift only).

- Simplify or remove the StokesParameter enum. While clearer, a simple string is likely sufficient.
