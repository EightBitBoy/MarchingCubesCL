#define STRINGIFY(A) #A // awesome hack to create a string from source code

std::string getKernelSource()
{
	std::string s = STRINGIFY(
		__kernel void marchingCubes(__global int *OUT)
		{
			int i = get_global_id(0);

			OUT[i] = i;
		}
	);

	return s;
}
