#include <stdexcept>
#include <utility>
#include <fstream>
#include <boost/thread/mutex.hpp>
#include <fantom/algorithm.hpp>
#include <fantom/graphics.hpp>
#include <fantom/fields.hpp>
#include <fantom/gui.hpp>
#include <CL/cl.hpp>
#include "MarchingCubesKernel.hpp"
#include "MarchingCubesTables.hpp"

#define __CL_ENABLE_EXCEPTIONS

using namespace std;
using namespace boost;
using namespace fantom;

namespace MC
{
	class MarchingCubes: public Algorithm
	{
	public:
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

		void printError(cl_int error, string description)
		{
			if(error != CL_SUCCESS)
			{
				debugLog() << description << " ERROR: " << error << endl;
			}
		}

		struct OptionsWindow
		{

			DockWindow mWindow;
			BoxLayout mLayout;
			MarchingCubes& mAlgo;
			Slider mSlider;
			

			OptionsWindow(MainWindow& mainWindow, MarchingCubes& algo)
				: mWindow(mainWindow, DockWindow::FFREE, "algorithm window"),
				mLayout(mWindow.getWidgetHolder(), false),
				mSlider(mLayout.addWidgetHolder(), 10, true, bind(&OptionsWindow::startPolygonizing, this)),
				mAlgo(algo)
			{
			}
			
			void startPolygonizing()
			{
				mAlgo.scheduleJob(bind(&MarchingCubes::polygonizeOpenCL, &mAlgo, mSlider.get()));
			}
		};

		void polygonizeOpenCL(size_t& value)
		{
			float isoValue = value;
			unique_lock<mutex> lock(mMutex);
			infoLog() << isoValue << endl;
		}

		void polygonizeNoOpenCL()
		{
		}

		mutex mMutex;
		Window<OptionsWindow> mWindow;
		unique_ptr<Graphics> polygonGroup;

		float* values;
		cl_float4* pointsVec;

		MarchingCubes(const Parameters& parameters): Algorithm(parameters), mWindow(*this)
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
			//polygonGroup->primitive().addSphere(Point3(0, 0, 0), 0.5, color);

			size_t numCells = grid->numCells();
			//size_t numCells = 1000;
			const size_t numCellPoints = 8;
			size_t numValues = numCells * numCellPoints;
			int time = 0;

			// load all scalar values into an array
			values = new float[numValues];
			for(Progress i(*this, "load values", numCells); i < numCells; ++i)
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
			debugLog() << "loading values finished" << endl;

			// load all points into an array
			pointsVec = new cl_float4[numValues];
			for(Progress i(*this, "load cells", numCells); i < numCells; ++i)
			{
				Cell cell = grid->cell(i);
				for(size_t j = 0; j < numCellPoints; ++j)
				{
					Point3 p = points[cell.index(j)];
					cl_float4 pVec;

					pVec.s[0] = p[0];
					pVec.s[1] = p[1];
					pVec.s[2] = p[2];

					pointsVec[(i * numCellPoints) + j] = pVec;
				}
			}
			debugLog() << "loading points finished" << endl;

			// ==================
			// CPU implementation
			// ==================
			if(useOpenCL == false)
			{
				debugLog() << ">>>> not using OpenCL <<<<" << endl;

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
				delete[] indices;
			}

			// =====================
			// OpenCL implementation
			// =====================
			if(useOpenCL == true)
			{
				debugLog() << ">>>> using OpenCL <<<<" << endl;

				cl_int error = CL_SUCCESS;
				cl_ulong infoNumber = 0;
				string infoString = "";

				// get the platforms
				vector<cl::Platform> platforms;
				cl::Platform::get(&platforms);

				for(int i = 0; i < platforms.size(); i++)
				{
					platforms[i].getInfo(CL_PLATFORM_NAME, &infoString);
					debugLog() << "platform #" << i << " name: " << infoString << endl;
					platforms[i].getInfo(CL_PLATFORM_VERSION, &infoString);
					debugLog() << "platform #" << i << " version: " << infoString << endl;
				}

				// create a context
				cl_context_properties properties[3] = { 
					CL_CONTEXT_PLATFORM, 
					(cl_context_properties)(platforms[0])(), 
					0 
				};
				cl::Context context(CL_DEVICE_TYPE_GPU, properties, NULL, NULL, &error);
				printError(error, "context");

				// get the devices
				vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
				if(devices.size() == 0)
				{
					throw runtime_error("No devices found!");
				}
				for(int i = 0; i < devices.size(); i++)
				{
					devices[i].getInfo(CL_DEVICE_NAME, &infoString);
					debugLog() << "device #"<< i << " name: " << infoString << endl;
					devices[i].getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &infoNumber);
					debugLog() << "device #"<< i << " global memory size (MB): " << (infoNumber/1024/1024) << endl;
					devices[i].getInfo(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, &infoNumber);
					debugLog() << "device #"<< i << " global memory cache size (MB): " << (infoNumber/1024/1024) << endl;
					devices[i].getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &infoNumber);
					debugLog() << "device #"<< i << " max clock frequency: " << infoNumber << endl;
					devices[i].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &infoNumber);
					debugLog() << "device #"<< i << " max compute units: " << infoNumber << endl;
					devices[i].getInfo(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, &infoNumber);
					debugLog() << "device #"<< i << " max constant buffer size (KB): " << (infoNumber/1024) << endl;
				}

				cl::CommandQueue queue = cl::CommandQueue(context, devices[0]);

				// prepare the kernel
				string s = getKernelSource();
				cl::Program::Sources source(1, std::make_pair(s.c_str(), s.length()+1));

				cl::Program program = cl::Program(context, source);
				error = program.build(devices);
				printError(error, "build");

				program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &infoString);
				debugLog() << "build info: " << infoString << endl;

				cl::Kernel kernel(program, "marchingCubes");

				// prepare buffers
				int* triTableArray = new int[256 * 16];
				for(int i = 0; i < 256; ++i)
				{
					for(int j = 0; j < 16; ++j)
					{
						triTableArray[(i * 16) + j] = triTable[i][j];
					}
				}

				cl::Buffer bufferEdgeTable = cl::Buffer(context, CL_MEM_READ_ONLY, 256 * sizeof(int));
				error = queue.enqueueWriteBuffer(bufferEdgeTable, CL_TRUE, 0, 256 * sizeof(int), edgeTable);
				printError(error, "bufferEdgeTable");

				cl::Buffer bufferTriTable = cl::Buffer(context, CL_MEM_READ_ONLY, 256 * 16 * sizeof(int));
				error = queue.enqueueWriteBuffer(bufferTriTable, CL_TRUE, 0, 256 * 16 * sizeof(int), triTableArray);
				printError(error, "bufferTriTable");

				cl::Buffer bufferValues = cl::Buffer(context, CL_MEM_READ_ONLY, numValues * sizeof(float));
				error = queue.enqueueWriteBuffer(bufferValues, CL_TRUE, 0, numValues * sizeof(float), values);
				printError(error, "bufferValues");

				cl::Buffer bufferPointsVec = cl::Buffer(context, CL_MEM_READ_ONLY, numValues *sizeof(cl_float4));
				error = queue.enqueueWriteBuffer(bufferPointsVec, CL_TRUE, 0, numValues *sizeof(cl_float4), pointsVec);
				printError(error, "bufferPointsVec");

				cl::Buffer bufferTriPoints = cl::Buffer(context, CL_MEM_WRITE_ONLY, 12 * numCells * sizeof(cl_float4));

				cl::Buffer bufferIndices = cl::Buffer(context, CL_MEM_WRITE_ONLY, numCells * sizeof(int));

				cl::Buffer bufferFloatTest = cl::Buffer(context, CL_MEM_WRITE_ONLY, numCells * sizeof(float));
				cl::Buffer bufferIntTest = cl::Buffer(context, CL_MEM_WRITE_ONLY, numCells * sizeof(int));


				kernel.setArg(0, isoValue);
				kernel.setArg(1, bufferEdgeTable);
				kernel.setArg(2, bufferTriTable);
				kernel.setArg(3, bufferValues);
				kernel.setArg(4, bufferPointsVec);
				kernel.setArg(5, bufferTriPoints);
				kernel.setArg(6, bufferIndices);

				kernel.setArg(7, bufferFloatTest);
				kernel.setArg(8, bufferIntTest);

				// run the kernel
				cl::NDRange global(numCells);
				cl::NDRange local(1);
				queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);

				// get the results
				float* floatTest = new float[numCells];
				queue.enqueueReadBuffer(bufferFloatTest, CL_TRUE, 0, numCells * sizeof(float), floatTest);
				for(int i = 0; i < numCells; i++)
				{
					//debugLog() << floatTest[i] << endl;
				}

				int* intTest = new int[numCells];
				queue.enqueueReadBuffer(bufferIntTest, CL_TRUE, 0, numCells * sizeof(int), intTest);
				for(int i = 0; i < numCells; i++)
				{
					//debugLog() << intTest[i] << endl;
				}
				
				cl_float4* triPoints = new cl_float4[12 * numCells];
				queue.enqueueReadBuffer(bufferTriPoints, CL_TRUE, 0, 12 * numCells * sizeof(cl_float4), triPoints);

				int* indices = new int[numCells];
				queue.enqueueReadBuffer(bufferIndices, CL_TRUE, 0, numCells * sizeof(int), indices);


				// draw the triangles
				for(size_t i = 0; i < numCells; ++i)
				{
					int index = indices[i];
					for(int j = 0; triTable[index][j] != -1; j = j+3)
					{
						vector<cl_float4> triPointsVec;
						triPointsVec.push_back(triPoints[(12 * i) + (triTable[index][j + 0])]);
						triPointsVec.push_back(triPoints[(12 * i) + (triTable[index][j + 1])]);
						triPointsVec.push_back(triPoints[(12 * i) + (triTable[index][j + 2])]);

						vector<Point3> triPoints;
						for(int k = 0; k < 3; ++k)
						{
							Point3 p;
							p[0] = triPointsVec[k].s[0];
							p[1] = triPointsVec[k].s[1];
							p[2] = triPointsVec[k].s[2];
							triPoints.push_back(p);
						}

						polygonGroup->primitive().add(Primitive::TRIANGLES).setColor(color).setVertices(triPoints);
					}
				}

				// free memory
				delete[] floatTest;
				delete[] intTest;
				delete[] triPoints;
				delete[] indices;
			}

			// free memory
			delete[] values;
			delete[] pointsVec;
		}
	};

	AlgorithmRegister<MarchingCubes> dummy("MarchingCubesCL", "executes MarchingCubesCL");
}
