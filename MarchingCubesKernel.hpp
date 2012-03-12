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
			__global float4* pointsVec,
			__global float4* triPoints,
			__global int* triNumber,

			__global float* floatTest,
			__global int* intTest
		)
		{
			int i = get_global_id(0);

			int index = 0;
			if (values[(i * 8) + 0] < isoValue) index |=   1;
			if (values[(i * 8) + 1] < isoValue) index |=   2;
			if (values[(i * 8) + 2] < isoValue) index |=   4;
			if (values[(i * 8) + 3] < isoValue) index |=   8;
			if (values[(i * 8) + 4] < isoValue) index |=  16;
			if (values[(i * 8) + 5] < isoValue) index |=  32;
			if (values[(i * 8) + 6] < isoValue) index |=  64;
			if (values[(i * 8) + 7] < isoValue) index |= 128;






			floatTest[i] = pointsVec[0].x;
			intTest[i] = index;
		}
	);

	return s;
}
