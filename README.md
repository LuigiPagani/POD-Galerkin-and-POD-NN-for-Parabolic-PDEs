# Reduced Order Models for Parabolic PDEs using POD-Galerkin and POD-NN Methods

Welcome to this repository, where we explore the development of reduced-order models for solving parabolic partial differential equations (PDEs). The work focuses on two distinct methods: the **POD-Galerkin** approach and the **POD-NN** (Proper Orthogonal Decomposition with Neural Networks) approach.

## What's Inside

This repository contains two Jupyter notebooks, each dedicated to one of the methods:

1. **POD-Galerkin Method**: This notebook implements the classical POD-Galerkin method, where the solution of the PDE is approximated using a reduced basis obtained from the Proper Orthogonal Decomposition. The Galerkin projection is then used to derive a system of ordinary differential equations (ODEs) that approximate the original PDE.

2. **POD-NN Method**: In this notebook, the POD method is combined with a neural network approach. The idea is to use a neural network to learn the dynamics of the reduced-order model, which can potentially offer better performance in terms of accuracy and computational efficiency compared to the traditional POD-Galerkin method.

## Getting Started

To explore these methods, simply open the notebooks in your Jupyter environment and follow along with the code and explanations. The notebooks are designed to be self-contained, with all the necessary steps to replicate the results and gain a deeper understanding of each method.
