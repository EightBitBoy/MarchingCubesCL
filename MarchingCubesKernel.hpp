#define STRINGIFY(A) #A // awesome hack to create a string from source code

std::string getKernelSource()
{
	std::string s = STRINGIFY(
		float4 interpolate(
		)
		{
		}

		__kernel void marchingCubes(
			float isoValue,
			__global int* edgeTable,
			__global int* triTable,
			__global float* values,

			__global float* floatTest,
			__global int* intTest
		)
		{
			int i = get_global_id(0);



			floatTest[i] = values[1];
			intTest[i] = triTable[96];
		}
	);

	return s;
}
