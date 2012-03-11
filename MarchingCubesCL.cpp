#include <stdexcept>
#include <utility>
#include <fstream>
#include <fantom/algorithm.hpp>
#include <fantom/graphics.hpp>
#include <fantom/fields.hpp>
#include <fantom/gui.hpp>
#include <CL/cl.hpp>
#include "MarchingCubesKernel.hpp"
#include "MarchingCubesTables.hpp"

#define __CL_ENABLE_EXCEPTIONS

using namespace std;
using namespace fantom;

namespace MC
{
	class MarchingCubes: public Algorithm
	{
	public:
		unique_ptr<Graphics> polygonGroup;

		struct Options: public Algorithm::Options
		{
			Options()
			{
				add<const TensorField<3, Scalar>*>("field", "", 0);
				add<const Grid<3>*>("grid", "", 0);
				add<bool>("use OpenCL", "switch between OpenCL (GPU) and CPU implementations", false);
				add<float>("iso value", "", 1.0);
				add<Color>("color", "", Color(0.75, 0.0, 0.0));
			}
		};

		Point3 interpolate(float isoValue, Point3 pointA, Point3 pointB, float valueA, float valueB)
		{
			float factor;
			Point3 result;

			factor = (isoValue - valueA) / (valueB - valueA);
			result[0] = pointA[0] + factor * (pointB[0] - pointA[0]);
			result[1] = pointA[1] + factor * (pointB[1] - pointA[1]);
			result[2] = pointA[2] + factor * (pointB[2] - pointA[2]);

			return result;
		}

		MarchingCubes(const Parameters& parameters): Algorithm(parameters)
		{
			// get the data
			const TensorField<3, Scalar>* field = parameters.get<const TensorField<3, Scalar>*>("field");
			const Grid<3>* grid = parameters.get<const Grid<3>*>("grid");
			const bool useOpenCL = parameters.get<bool>("use OpenCL");
			const float isoValue = parameters.get<float>("iso value");
			const Color color = parameters.get<Color>("color");

			if(field == false)
			{
				throw runtime_error("field not set!");
			}
			if(grid == false)
			{
				throw runtime_error("grid not set!");
			}
 
			auto evaluator = field->makeEvaluator();
			auto& points = grid->parent().points();

			polygonGroup = makeGraphics("iso surface");
			polygonGroup->primitive().addSphere(Point3(0, 0, 0), 0.5, color);
			// ==================
			// CPU implementation
			// ==================
			if(useOpenCL == false)
			{
				// TODO use this!
				//size_t numCells = grid->numCells();
				size_t numCells = 1000;
				const size_t numCellPoints = 8;
				size_t numValues = numCells * numCellPoints;
				int time = 0;

				// load all scalar values into an array
				float* values = new float[numValues];
				for(Progress i(*this, "load data", numCells); i < numCells; ++i)
				{
					Cell cell = grid->cell(i);
					for(size_t j = 0; j < numCellPoints; ++j)
					{
						// TODO USE THE CORRECT VALUES
						Point3 p = points[cell.index(j)];
						if(evaluator->reset(p, time))
						{
						}
						else
						{
							
						}
						values[(i * numCellPoints) + j] = (float)(rand() % 10);
					}
				}
				debugLog() << "loading finished" << endl;

				// get the cube indices
				int* indices = new int[numCells];
				for(Progress i(*this, "calculating indices", numCells); i < numCells; ++i)
				{
					int index = 0;

					if (values[(i * numCellPoints) + 0] < isoValue) index |=   1;
					if (values[(i * numCellPoints) + 1] < isoValue) index |=   2;
					if (values[(i * numCellPoints) + 2] < isoValue) index |=   4;
					if (values[(i * numCellPoints) + 3] < isoValue) index |=   8;
					if (values[(i * numCellPoints) + 4] < isoValue) index |=  16;
					if (values[(i * numCellPoints) + 5] < isoValue) index |=  32;
					if (values[(i * numCellPoints) + 6] < isoValue) index |=  64;
					if (values[(i * numCellPoints) + 7] < isoValue) index |= 128;

					indices[i] = index;
				}
				debugLog() << "indexing finished" << endl;


				// calculate and draw polygons
				for(Progress i(*this, "polygonizing", numCells); i < numCells; ++i)
				{
					Cell cell = grid->cell(i);
					Point3 cellPoints[numCellPoints];
					float cellValues[numCellPoints];

					for(int j = 0; j < numCellPoints; ++j)
					{
						cellPoints[j] = points[cell.index(j)];
						cellValues[j] = values[(i * numCellPoints) +j];
					}

					int index = indices[i];
					int indexValue = edgeTable[index];
					Point3 vertices[12];

					if(indexValue == 0)
						continue;
					if(indexValue &    1)
						vertices[ 0] = interpolate(isoValue, cellPoints[0], cellPoints[1], cellValues[0], cellValues[1]);
					if(indexValue &    2)
						vertices[ 1] = interpolate(isoValue, cellPoints[1], cellPoints[2], cellValues[1], cellValues[2]);
					if(indexValue &    4)
						vertices[ 2] = interpolate(isoValue, cellPoints[2], cellPoints[3], cellValues[2], cellValues[3]);
					if(indexValue &    8)
						vertices[ 3] = interpolate(isoValue, cellPoints[3], cellPoints[0], cellValues[3], cellValues[0]);
					if(indexValue &   16)
						vertices[ 4] = interpolate(isoValue, cellPoints[4], cellPoints[5], cellValues[4], cellValues[5]);
					if(indexValue &   32)
						vertices[ 5] = interpolate(isoValue, cellPoints[5], cellPoints[6], cellValues[5], cellValues[6]);
					if(indexValue &   64)
						vertices[ 6] = interpolate(isoValue, cellPoints[6], cellPoints[7], cellValues[6], cellValues[7]);
					if(indexValue &  128)
						vertices[ 7] = interpolate(isoValue, cellPoints[7], cellPoints[4], cellValues[7], cellValues[4]);
					if(indexValue &  256)
						vertices[ 8] = interpolate(isoValue, cellPoints[0], cellPoints[4], cellValues[0], cellValues[4]);
					if(indexValue &  512)
						vertices[ 9] = interpolate(isoValue, cellPoints[1], cellPoints[5], cellValues[1], cellValues[5]);
					if(indexValue & 1024)
						vertices[10] = interpolate(isoValue, cellPoints[2], cellPoints[6], cellValues[2], cellValues[6]);
					if(indexValue & 2048)
						vertices[11] = interpolate(isoValue, cellPoints[3], cellPoints[7], cellValues[3], cellValues[7]);


					for(int j = 0; triTable[index][j] != -1; j = j+3)
					{
						vector<Point3> triPoints;
						triPoints.push_back(vertices[triTable[index][j    ]]);
						triPoints.push_back(vertices[triTable[index][j + 1]]);
						triPoints.push_back(vertices[triTable[index][j + 2]]);

						polygonGroup->primitive().add(Primitive::TRIANGLES).setColor(color).setVertices(triPoints);
					}
				}
				debugLog() << "polygonization finished" << endl;

				// free memory
				delete[] values;
				delete[] indices;
			}

			// =====================
			// OpenCL implementation
			// =====================
			if(useOpenCL == true)
			{
				throw runtime_error("OpenCL algorithm not yet implemented!");
			}
		}
	};

	AlgorithmRegister<MarchingCubes> dummy("MarchingCubesCL", "executes MarchingCubesCL");
}
