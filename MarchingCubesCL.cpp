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
			}
		};

		Point3 interpolate(float isoValue, Point3 pointA, Point3 pointB, float valueA, float valueB)
		{
			float factor;
			Point3 result;
		}

		MarchingCubes(const Parameters& parameters): Algorithm(parameters)
		{
			// get the data
			const TensorField<3, Scalar>* field = parameters.get<const TensorField<3, Scalar>*>("field");
			const Grid<3>* grid = parameters.get<const Grid<3>*>("grid");
			const bool useOpenCL = parameters.get<bool>("use OpenCL");
			const float isoValue = parameters.get<float>("iso value");

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

			// ==================
			// CPU implementation
			// ==================
			if(useOpenCL == false)
			{
				// TODO use this!
				//size_t numCells = grid->numCells();
				size_t numCells = 6000;
				size_t numCellPoints = 8;
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
					debugLog() << index << endl;
				}
				debugLog() << "indexing finished" << endl;

				// calculate and draw polygons
				for(Progress i(*this, "polygonizing", numCells); i < numCells; ++i)
				{
					int indexValue = edgeTable[indices[i]];
					if(indexValue == 0)
						continue;
					
					//if(indexValue & 1)
					//	interpolate();
				}
				debugLog() << "polygonizing finished" << endl;

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
