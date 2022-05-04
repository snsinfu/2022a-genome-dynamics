Model features:

- Anaphase initialization at 10Mb resolution, plus spherical packing
- Spline refinement of 10Mb conformation to 100kb one
- Interphase simulation with:
  - A/B/u chromatin beads
  - Nucleolus formation
  - Indenpendently moving nuclear semiaxes under external compression


### Trajectory file format

HDF5 with the following hierarchy:

```
/
  metadata/
    config              string
    particle_types      int (N)
    ab_factors          float (N, 2)
    chromosome_ranges   int (*, 2)
    centromere_ranges   int (*, 2)
    nucleolus_ranges    int (*, 2)
    nucleolus_bonds     int (*, 2)

  snapshots/
    spindle/
      0/
        positions           float (Ni, 3)
      ...

    packing/
      0/
        positions           float (Ni, 3)
      ...

    relaxation/
      0/
        context             string
        positions           float (N, 3)
      ...

    interphase/
      0/
        context             string
        positions           float (N, 3)
        contact_map         int (*, 3)
      ...
```
