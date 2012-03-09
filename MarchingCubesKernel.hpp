#define STRINGIFY(A) #A // awesome hack to create a string from source code

std::string getKernelSource()
{
	std::string s = STRINGIFY(
		__kernel void vector_add(__global const int *A, const int *B, __global int *C, __global const int ISO, __constant int *EDGES, __global const float *VALUES, __global float *OUTPUT, __global int *INDICES) {
			int i = get_global_id(0);
			
			int cubeIndex = 0;
			if(VALUES[(8 * i) + 0] < ISO) cubeIndex |= 1;
			if(VALUES[(8 * i) + 1] < ISO) cubeIndex |= 2;
			if(VALUES[(8 * i) + 2] < ISO) cubeIndex |= 4;
			if(VALUES[(8 * i) + 3] < ISO) cubeIndex |= 8;
			if(VALUES[(8 * i) + 4] < ISO) cubeIndex |= 16;
			if(VALUES[(8 * i) + 5] < ISO) cubeIndex |= 32;
			if(VALUES[(8 * i) + 6] < ISO) cubeIndex |= 64;
			if(VALUES[(8 * i) + 7] < ISO) cubeIndex |= 128;
			INDICES[i] = cubeIndex;
			

			OUTPUT[i] = ISO;
			if(i == 5)
			{
				C[i] = 999;
			}
			else
			{
				C[i] = A[i] + B[i];
			}
		}
	);

	return s;
}
