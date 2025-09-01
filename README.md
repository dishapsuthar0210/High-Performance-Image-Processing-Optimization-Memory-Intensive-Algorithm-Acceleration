# High-Performance-Image-Processing-Optimization-Memory-Intensive-Algorithm-Acceleration
---
## Setup and Execution

Follow these steps to set up and run the project:

1. **Clone the repository** into your local machine or download the `.tar` file and extract it into your working directory.  
2. **Windows users**: Install [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows Subsystem for Linux) to run a virtual Linux environment and create an account.  
3. In WSL, **navigate** to the directory where you placed the project files.  
4. The source code will be linked with the object code we supply into a driver binary. To create this binary, execute the following commands:

   ```bash
   unix> make driver
---
5. To test your implementations, run:
```bash
unix> ./driver
```
- Achieved 26.4x speedup for image rotation algorithms through advanced optimization techniques including cache-friendly blocking (32x32 tiles), loop unrolling, pointer arithmetic optimization, and memory access pattern restructuring
- Delivered 68.3x performance improvement for image smoothing operations by implementing case-split boundary handling, eliminating redundant computations, and applying sliding window algorithms for 3x3 convolution kernels
- Optimized memory-intensive C code for real-time image processing operations (90° rotation and Gaussian blur) on square matrices up to 1024x1024 pixels, reducing cycles-per-element (CPE) from 94.5 to 4.4 for rotation and 722 to 10.5 for smoothing
- Applied low-level optimization strategies including cache locality improvements, computational complexity reduction, strength reduction techniques, and eliminated function call overhead through strategic code inlining
