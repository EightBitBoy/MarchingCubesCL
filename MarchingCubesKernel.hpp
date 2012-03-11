#define STRINGIFY(A) #A // awesome hack to create a string from source code

std::string getKernelSource()
{
	std::string s = STRINGIFY(
		__kernel void marchingCubes(
			float isoValue,
			__global float* values,
			__global float4* pointsVec,
			__global float* OUTTEST
		)
		{
			int i = get_global_id(0);

			OUTTEST[i] = pointsVec[i].x;
		}
	);

	return s;
}
