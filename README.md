# SecONNds: Secure Outsourced Neural Network Inference on ImageNet

## Introduction

SecONNds is a framework for secure outsourced neural network inference, particularly focusing on ImageNet models. It allows a client to perform inference on a server's model without revealing the client's input or the server's model weights.

## Prerequisites

Ensure that the following tools are installed on your system:

- `g++`
- `make`
- `git`
- `cmake` (version 3.10 or higher)

## Building the Project

To build SecONNds, run:

```bash
bash scripts/build.sh [options]
```

### Build Options

- `-clean`: Removes the build directory before building. **Default:** Do not clean.
- `-debug`: Builds in Debug mode. **Default:** Release mode.
- `--trip_trials`: Runs triple generation trials. **Default:** Off.
- `--track_he_noise` or `-noise`: Tracks Homomorphic Encryption (HE) noise during computation. **Default:** Off.
- `--verify_layerwise` or `-verify`: Verifies layer-wise outputs during execution. **Default:** Off.

### Examples

- **Clean build in Release mode:**

  ```bash
  bash scripts/build.sh -clean
  ```

- **Build in Debug mode with HE noise tracking:**

  ```bash
  bash scripts/build.sh -debug -noise
  ```

### Dependencies

The build script will automatically check for and build necessary dependencies if they are not present. These include:

- [Eigen](https://gitlab.com/libeigen/eigen)
- [EMP-Toolkit](https://github.com/emp-toolkit) (`emp-tool`, `emp-ot`)
- [Zstandard](https://github.com/facebook/zstd) (`zstd`)
- [Intel HEXL](https://github.com/intel/hexl) (`hexl`)
- [Microsoft SEAL](https://github.com/microsoft/SEAL) (`SEAL`)
- [SEAL-CUDA](https://github.com/tnishimoto/SEAL-Embedded) (`seal-cuda`)

## Running the Project

To run SecONNds, use:

```bash
bash scripts/run-cnn.sh [server|client] [framework] [network] [options]
```

### Required Arguments

1. **Role**: `server` or `client`
2. **Framework**: Choose from:
   - `seconnds_2`
   - `seconnds_p`
   - `cheetah`
   - `SCI_HE`
3. **Network**: Choose from:
   - `sqnet`
   - `resnet50`
   - `densenet121`

### Optional Arguments

- `--mill_low_rnd` or `-mlr`: Use Millionaires' protocol with lower rounds. **Default:** Off.
  - *Note:* Applicable for `seconnds_2` and `seconnds_p` frameworks.
- `--conv_ntt` or `-ntt`: Use convolution with NTT preprocessing. **Default:** Off.
  - *Note:* Applicable for `cheetah` and `SCI_HE` frameworks.
- `-bl=[bit length]`: Set the secret sharing bit length. **Default:** `32`.
  - *Note:* Must be the same for both server and client.
- `-j=[number of threads]`: Set the number of threads. **Default:** `1`.
  - *Note:* Must be the same for both server and client.
- `-debug`: Run in Debug mode using GDB. **Default:** Off.
- `--no_log` or `-nl`: Do not generate a log file. **Default:** Log files are generated.
- `-l=[log file number]`: Prefix the log file name with a number. **Default:** No prefix.

### Examples

- **Run as server using `seconnds_2` on `sqnet` network:**

  ```bash
  bash scripts/run-cnn.sh server seconnds_2 sqnet
  ```

- **Run as client using `cheetah` on `resnet50` with NTT preprocessing:**

  ```bash
  bash scripts/run-cnn.sh client cheetah resnet50 -ntt
  ```

- **Run as server with 4 threads and secret sharing bit length 37:**

  ```bash
  bash scripts/run-cnn.sh server seconnds_p sqnet -j=4 -bl=37
  ```

### Notes on Options

- **`-bl` Option:**
  - Sets the bit length for secret sharing.
  - Must be consistent between server and client.
- **`-j` Option:**
  - Sets the number of threads for execution.
  - Must be consistent between server and client.
- **`-mlr` Option:**
  - Activates the Millionaires' protocol with lower rounds.
  - Applicable only to `seconnds_2` and `seconnds_p` frameworks.
  - Must be consistent between server and client.
- **`-ntt` Option:**
  - Enables convolution with NTT preprocessing.
  - Applicable only to `cheetah` and `SCI_HE` frameworks.

## Running All Frameworks

To run all frameworks for a given network, use:

```bash
bash scripts/run-cnn-all-fw.sh [server|client] [network] [options]
```

- This script sequentially runs the specified network with all frameworks (`seconnds_2`, `seconnds_p`, `cheetah`, `SCI_HE`) and corresponding options.

### Example

- **Run all frameworks as client on `densenet121` network:**

  ```bash
  bash scripts/run-cnn-all-fw.sh client densenet121
  ```

## Common Configuration

The scripts use common configurations defined in `scripts/common.sh`. Important variables include:

- **`SERVER_IP`**: IP address of the server.
  - **Default:** `127.0.0.1`
- **`SERVER_PORT`**: Port number for communication.
  - **Default:** `12345`
- **`FXP_SCALE`**: Fixed-point scaling factor.
  - **Default:** `12`
- **`SS_BITLEN`**: Secret sharing bit length.
  - **Default:** `32`
- **`NUM_THREADS`**: Number of threads for execution.
  - **Default:** `1`
- **`CSIZE`**: Chunk size for triple generation.
  - **Default:** `32000`

*Ensure that `SS_BITLEN` and `NUM_THREADS` are the same on both server and client.*

## Logs

- By default, logs are saved in the `logs/` directory.
- The log file name format is:

  ```
  [network]-bl[bit length]-j[threads]-[framework]-[options]-[role].log
  ```

- Use `-l=[number]` to prefix the log file with a number.

### Disable Logging

- Use the `--no_log` or `-nl` option to prevent log file generation.

## Debugging

- Use the `-debug` option to run the program under GDB for debugging purposes.

  ```bash
  bash scripts/run-cnn.sh [server|client] [framework] [network] -debug
  ```

## Additional Information

- **Dependencies are patched and built automatically.**
  - Custom patches are applied to some dependencies to ensure compatibility.
- **Ensure consistent options between server and client.**
  - Options like `-mlr`, `-bl`, and `-j` must be the same on both sides.
- **Data Directory:**
  - Input files and model weights are expected to be in the `pretrained/` directory.
  - The scripts will create a `data/` directory to store intermediate data if needed.

## Troubleshooting

- **CMake Version:**
  - If you encounter issues during the build, ensure that your CMake version is 3.10 or higher.
- **Missing Dependencies:**
  - If a dependency is missing, the build script should automatically attempt to build it.
  - If problems persist, manually check the `deps/` and `build/` directories.

## Contact

For issues or questions, please contact the authors of the repository.
