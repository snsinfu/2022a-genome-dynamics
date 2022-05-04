analyze_lamina distance <output> <inputs>...
  Takes simulation trajectories and outputs wall-distance information to
  `/distance/{input}` and `/normalized_distance/{input}`.

analyze_lamina contact <output>
  Takes analysis file, reads `/normalized_distance` datasets and outputs contact
  information to `/contact/{method}/{input}` and `/average_contact/{method}`.

  options:
    --contact-model (naive|tanh|gaussian)
    --tanh-delta 0.1
    --gaussian-sigma 2.0
