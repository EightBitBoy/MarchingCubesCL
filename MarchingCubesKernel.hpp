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
			__global float* floatTest,
			__global int* intTest
		)
		{
			int i = get_global_id(0);



			floatTest[i] = isoValue;
			intTest[i] = i;
		}
	);

	return s;
}
