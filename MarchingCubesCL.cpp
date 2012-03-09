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
			}
		};

		Point3 interpolate(Point3 pointA, Point3 pointB, float valueA, float valueB, float isoValue)
		{
			Point3 result;

			float factor;
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

				size_t numCells = grid->numCells();
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
						Point3 p = points[cell.index(j)];
						if(evaluator->reset(p, time))
						{
						}
						else
						{
						}
					}
				}

				// get the cube indices
				int* indices = new int[numCells];
				for(Progress i(*this, "calculating indices", numCells); i < numCells; ++i)
				{
					for(size_t j = 0; j < numCellPoints; ++j)
					{
						
					}
				}

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
