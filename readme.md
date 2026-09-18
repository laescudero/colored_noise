# Zeros of the Spectrogram of Colored Noise

This repository contains the numerical implementation and experiments described in the paper **"Zeros of the Spectrogram of Colored Noise"** by L. A. Escudero, G. Koliander, and J. L. Romero.

The code provides tools to simulate the Short-Time Fourier Transform (STFT) of signals embedded in complex Gaussian colored noise, extract its zeros, and compute the estimators described in the manuscript. 

## Repository Structure

The repository is organized around two main executable scripts:

* **`main.m`**: The core experimental pipeline. It runs the simulations from our paper.
* **`demo_spectrogram.m`**: A standalone script to generate the spectrogram and the distribution of its zeros for a deterministic signal embedded in colored noise. 
	
## Requirements and Compatibility

The code has been designed with cross-compatibility in mind. You need either MATLAB or Octave.

* **MATLAB**: The Parallel Computing Toolbox is highly recommended to speed up the simulations.
* **GNU Octave**: For Parallelization we recommend installing the `parallel` package.

## How to Run

1. Clone this repository to your computer.
2. Open your MATLAB or Octave environment and set the repository folder as your working directory.
3. For the experiments described in Section 5 and 6 of our paper, execute the main function.
4. Choose Option 1, 2 or 3 to choose the experiment.
5. Choose Option 1 to only compute the realizations. Option 2 to produce the plots, if you already have ran Option 1. Option 3 to compute and plot at once.
6. Alternatively, to visualize a single realization, simply run demo_spectrogram.

The produced plots are saved in **`images/`**. Alternatively, you can set **`export_raw_files = true;`** in main.m to obtain the raw files used in the plots in our article. Then, you have to compile the `*.tex` files included in the repository.

By default, the code will reproduce our results (if you run it on Octave 9.4.0). If you want to change this configuration, you have to set **`experimentsParam.use_fixed_seed = false`** in main.m.

## Citation
If you use this code in your research, please cite our paper:

```
@misc{escudero2026zeros,
  title        = {Zeros of the Spectrogram of Colored Noise},
  author       = {Escudero, Luis Alberto and Koliander, G{\"u}nther and Romero, Jos{\'e} Luis},
  year         = {2026},
  howpublished = {arXiv preprint, \url{https://arxiv.org/abs/2609.18689}}
}
```

## License
This project is open-source and available under the MIT License.
