# MEnDER-1D
Multi-Energy Differential Evolution Reconstruction (MEnDER 1D) for Proton Deflectometry

This code is designed to reconstruct magnetic field deflections, and thus the path-integrated magnetic field, from sets of proton images at two distinct proton probe energies. A differential evolution (DE) algorithm is used to iteratively update a population of solution candidates of the magnetic deflections of the protons for reconstructing the input images, selecting improved candidates as they are discovered. This algorithm was written using MATLAB (R2019a) and makes use of the Image Processing and Parallel Processing Toolboxes. This release contains 7 relevant files:

MEnDER_1D_Run.m is a script to set up and run the algorithm on a set of test problems.

Ellipsoid_Initialize.m, Wire_Initialize.m, and Truncated_Wire_Initialize.m, called from the run file, are functions used to generate proton images from specific test magnetic field structures.

MEnDER_1D_DE.m is the main function to reconstruct a set of input proton images using DE, calling the NewPopulation.m function to generate populations of solution candidates. 

radialweightedHist.m is a simple function to generate the final pseudo-1D proton intensity maps (weighted by initial and final image position) within MEnDER_1D_DE.m.


These tools were used to generate the results displayed in: J. M. Levesque and L. J. Beesley, “Reconstructing magnetic deflections from sets of proton images using differential evolution,” RSI 2021.


© 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 
This program is open source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
