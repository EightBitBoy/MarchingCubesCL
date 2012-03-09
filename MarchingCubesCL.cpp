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
				add<float>("isoValue", "the iso value used by the Marching Cubes algorithm", 0.00080);
				add<Color>("color", "", Color(0.75, 0.0, 0.0));
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
			const float isoValue = parameters.get<float>("isoValue");
			const Color color = parameters.get<Color>("color");

			if(field == false)
			{
				throw runtime_error("field not set!");
			}
			if(grid == false)
			{
				throw runtime_error("grid not set!");
			}

			// add some test objects
			polygonGroup = makeGraphics("polygonGroup");
			polygonGroup->primitive().addSphere(Point3(0, 0, 0), 0.5, color);
 
			try { 
				auto evaluator = field->makeEvaluator();
				auto& points = grid->parent().points();

				size_t numCells = grid->numCells();
				size_t numCellPoints = 8;
				size_t numValues = numCells * numCellPoints;
				int time = 0;

				// put all scalar values into an array
				float * values = new float[numValues];

				for(Progress i(*this, "loading data", numCells); i < numCells; ++i)
				{
					Cell cell = grid->cell(i);
					for(size_t j = 0; j < numCellPoints; ++j)
					{
						Point3 p = points[cell.index(j)];
						if(evaluator->reset(p, time))
						{
							values[(i * numCellPoints) + j] = evaluator->value()();
						}
						else
						{
							values[(i * numCellPoints) + j] = 0;
						}
					}
				}
				
				// free memory
				delete[] values;





				/*
				for(Progress i(*this, "working, size", numberOfCells); i < numberOfCells; ++i)
				{
					debugLog() << "working" << endl;
					Cell cell = grid->cell(i);

					Point3 cubePoints[8];
					for(size_t j = 0; j < numberOfCellPoints; ++j)
					{
						cubePoints[j] = points[cell.index(j)];
					}

					int index = cubeIndices[i];
					Point3 vertices[12];
					if(edgeTable[index] == 0)
						continue;
					

					if(edgeTable[index] & 1)
						vertices[0] = interpolate(cubePoints[0], cubePoints[1], values[(8 * i) + 0], values[(8 * i) + 1], isoValue);
					if(edgeTable[index] & 2)
						vertices[1] = interpolate(cubePoints[1], cubePoints[2], values[(8 * i) + 1], values[(8 * i) + 2], isoValue);
					if(edgeTable[index] & 4)
						vertices[2] = interpolate(cubePoints[2], cubePoints[3], values[(8 * i) + 2], values[(8 * i) + 3], isoValue);
					if(edgeTable[index] & 8)
						vertices[3] = interpolate(cubePoints[3], cubePoints[0], values[(8 * i) + 3], values[(8 * i) + 0], isoValue);
					if(edgeTable[index] & 16)
						vertices[4] = interpolate(cubePoints[4], cubePoints[5], values[(8 * i) + 4], values[(8 * i) + 5], isoValue);
					if(edgeTable[index] & 32)
						vertices[5] = interpolate(cubePoints[5], cubePoints[6], values[(8 * i) + 5], values[(8 * i) + 6], isoValue);
					if(edgeTable[index] & 64)
						vertices[6] = interpolate(cubePoints[6], cubePoints[7], values[(8 * i) + 6], values[(8 * i) + 7], isoValue);
					if(edgeTable[index] & 128)
						vertices[7] = interpolate(cubePoints[7], cubePoints[4], values[(8 * i) + 7], values[(8 * i) + 4], isoValue);
					if(edgeTable[index] & 256)
						vertices[8] = interpolate(cubePoints[0], cubePoints[4], values[(8 * i) + 0], values[(8 * i) + 4], isoValue);
					if(edgeTable[index] & 512)
						vertices[9] = interpolate(cubePoints[1], cubePoints[5], values[(8 * i) + 1], values[(8 * i) + 5], isoValue);
					if(edgeTable[index] & 1024)
						vertices[10] = interpolate(cubePoints[2], cubePoints[6], values[(8 * i) + 2], values[(8 * i) + 6], isoValue);
					if(edgeTable[index] & 2048)
						vertices[11] = interpolate(cubePoints[3], cubePoints[7], values[(8 * i) + 3], values[(8 * i) + 7], isoValue);

					for(int j = 0; triTable[index][j] != -1; j += 3)
					{
						debugLog() << "YEAH" << endl;
						vector<Point3> polygonPoints;
						polygonPoints.push_back(vertices[j + 0]);
						polygonPoints.push_back(vertices[j + 1]);
						polygonPoints.push_back(vertices[j + 2]);

						polygonGroup->primitive().add(Primitive::POLYGON).setColor(color).setVertices(polygonPoints);
					}
				}
				*/
			}
			catch(...)
			{
				throw runtime_error("SOMETHING WENT WRONG!!!");
			}
		}
	};

	AlgorithmRegister<MarchingCubes> dummy("MarchingCubesCL", "executes MarchingCubesCL");
}
