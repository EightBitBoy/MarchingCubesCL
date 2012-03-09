# MarchingCubesCL

## About
This project is about an implementation of the Marching Cubes algorithm using OpenCL. The result is an interactive plugin for the visualization system FAnToM.
More information at http://www.informatik.uni-leipzig.de/fantom/.

## Usage
- follow the instructions in the FAnToM manual to create a plugin project
- clone the repository and move all the files it contains to the project source folder, it should be possible to use the repository folder directly as source folder
- move the "CL" folder from the repository to the "include" folder
- move the file "libOpenCL.so" to the "lib" folder
- follow further manual instructions and run cmake, SET THE CORRECT INCLUDE AND INSTALL PATHS
- IMPORTANT: run cmake-gui and add to the variable "CMAKE_CXX_STANDARD_LIBRARIES" the argument "-lOpenCL" (without quotation marks of course)
- run make to compile the plugin