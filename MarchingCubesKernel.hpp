#define STRINGIFY(A) #A // awesome hack to create a string from source code

std::string getKernelSource()
{
	std::string s = STRINGIFY(
		float4 interpolate(
			float isoValue,
			float4 pointA,
			float4 pointB,
			float valueA,
			float valueB
		)
		{
			float factor;
			float4 result;

			factor = (isoValue - valueA) / (valueB - valueA);

			result.x = pointA.x + factor * (pointB.x - pointA.x);
			result.y = pointA.y + factor * (pointB.y - pointA.y);
			result.z = pointA.z + factor * (pointB.z - pointA.z);

			return result;
		}

		__kernel void marchingCubes(
			float isoValue,
			__global int* edgeTable,
			__global int* triTable,
			__global float* values,
			__global float4* pointsVec,
			__global float4* triPoints,
			__global int* indices
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
			indices[i] = index;

			if(edgeTable[index] &    1)
				triPoints[(i * 12) +  0] = interpolate(isoValue, pointsVec[(i * 8) + 0], pointsVec[(i * 8) + 1], values[(i * 8) + 0], values[(i * 8) + 1]);
			if(edgeTable[index] &    2)
				triPoints[(i * 12) +  1] = interpolate(isoValue, pointsVec[(i * 8) + 1], pointsVec[(i * 8) + 2], values[(i * 8) + 1], values[(i * 8) + 2]);
			if(edgeTable[index] &    4)
				triPoints[(i * 12) +  2] = interpolate(isoValue, pointsVec[(i * 8) + 2], pointsVec[(i * 8) + 3], values[(i * 8) + 2], values[(i * 8) + 3]);
			if(edgeTable[index] &    8)
				triPoints[(i * 12) +  3] = interpolate(isoValue, pointsVec[(i * 8) + 3], pointsVec[(i * 8) + 0], values[(i * 8) + 3], values[(i * 8) + 0]);
			if(edgeTable[index] &   16)
				triPoints[(i * 12) +  4] = interpolate(isoValue, pointsVec[(i * 8) + 4], pointsVec[(i * 8) + 5], values[(i * 8) + 4], values[(i * 8) + 5]);
			if(edgeTable[index] &   32)
				triPoints[(i * 12) +  5] = interpolate(isoValue, pointsVec[(i * 8) + 5], pointsVec[(i * 8) + 6], values[(i * 8) + 5], values[(i * 8) + 6]);
			if(edgeTable[index] &   64)
				triPoints[(i * 12) +  6] = interpolate(isoValue, pointsVec[(i * 8) + 6], pointsVec[(i * 8) + 7], values[(i * 8) + 6], values[(i * 8) + 7]);
			if(edgeTable[index] &  128)
				triPoints[(i * 12) +  7] = interpolate(isoValue, pointsVec[(i * 8) + 7], pointsVec[(i * 8) + 4], values[(i * 8) + 7], values[(i * 8) + 4]);
			if(edgeTable[index] &  265)
				triPoints[(i * 12) +  8] = interpolate(isoValue, pointsVec[(i * 8) + 0], pointsVec[(i * 8) + 4], values[(i * 8) + 0], values[(i * 8) + 4]);
			if(edgeTable[index] &  512)
				triPoints[(i * 12) +  9] = interpolate(isoValue, pointsVec[(i * 8) + 1], pointsVec[(i * 8) + 5], values[(i * 8) + 1], values[(i * 8) + 5]);
			if(edgeTable[index] & 1024)
				triPoints[(i * 12) + 10] = interpolate(isoValue, pointsVec[(i * 8) + 2], pointsVec[(i * 8) + 6], values[(i * 8) + 2], values[(i * 8) + 6]);
			if(edgeTable[index] & 2049)
				triPoints[(i * 12) + 11] = interpolate(isoValue, pointsVec[(i * 8) + 3], pointsVec[(i * 8) + 7], values[(i * 8) + 3], values[(i * 8) + 7]);
		}
	);

	return s;
}
