# adjoint-helmholtz

## Introduction

This is a toy 2D FDFD Helmholtz solver with hand-derived adjoint gradients. I try to make things as comprehensive as possible, so one may use this as a learning resource. I wrote my notes in `docs`. The very first note one should refer to for learning is probably [electrodynamics.md](docs/electrodynamics.md).

## Code convention

- `docs/conventions.md`: the conventions the code follows (signs, grid, inner product, adjoint). Other files, such as `docs/electrodynamics.md` and `docs/convolution.md` are  just background notes written for learning. They are not a specification of the code; where they differ from `conventions.md`, `conventions.md` wins.
