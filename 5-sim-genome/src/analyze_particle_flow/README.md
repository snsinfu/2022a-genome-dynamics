This module impelments a command for computing flow velocity in the simulated
system. There are two definitions of flow:

- Mean velocity of particles around each particle
- Mean velocity of particles in a bin around grid points


analyze_flow particle --velocity-delay 1 --scan-radius 0.6 output.h5 ...

analyze_flow grid --velocity-delay 1 --grid-interval 0.3 --scan-radius 0.6 output.h5 ...


```
/
  flow/
    particle_scan06/
      output-1  (701, 62178, 3)
      output-2  (701, 62178, 3)
      ...

    grid03_scan06/
      output-1  (701, *, *, *, 3)
      output-2  (701, *, *, *, 3)
      ...
```
