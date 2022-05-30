# Whole-genome chromatin dynamics simulation

Source code for the paper
**S. Fujishiro and M. Sasai "Generation of dynamic three-dimensional genome structure through phase separation of chromatin" (2022)**.
The project covers:

- Novel annotation of local chromatin using *NCI*, or the diagonal signal of Hi-C contact matrix
- 1kb-resolution chromatin dynamics simulations
- 100kb-resolution whole-genome interphase chromatin dynamics simulations
- 10Mb-resolution whole-genome ana/telophase chromatin dynamics simulations
- Coarse-graining of chromatin polymer based on the PRISM-theory

| The paper                   |
|-----------------------------|
| [Published version][pub]    |
| [bioRxiv preprint][bioRxiv] |

[bioRxiv]: https://www.biorxiv.org/content/10.1101/2021.05.06.443035
[pub]: https://doi.org/10.1073/pnas.2109838119

- [Directories](#directories)
- [License](#license)

## Directories

- [1-hic](1-hic): Preparation of chromosome contact maps. Downloads and preprocesses Hi-C dataset.
- [2-signal](2-signal): Analysis of the contact maps. Computes compartment signal and NCI.
- [3-sim-1kb](3-sim-1kb): Chromatin dynamics at 1kb resolution. Incorporates kinetic HP-1 attraction and cohesin looping. Estimates effective interaction among 100kb chromatin regions using PRISM theory.
- [4-sim-ab](4-sim-ab): Simulation of phase separation in purely repulsive polymer blend.
- [5-sim-genome](5-sim-genome): Simulation of human genome in ana/telophase and interphase.

## License

MIT License
