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
#define SENSITIVITY 1000

using namespace std;
using namespace boost;
using namespace fantom;

// TODO delete polygonGroup???

namespace MC
{
	class MarchingCubes: public Algorithm
	{
	public:
		struct Options: public Algorithm::Options
		{
			Options()
			{
				add<const TensorFieldGridBased<3, Scalar>*>("field", "", 0);
				add<const Grid<3>*>("grid", "", 0);
				add<float>("max iso value", "", 10.0);
				add<bool>("debug", "provides more output", false);
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
				mSlider(mLayout.addWidgetHolder(), SENSITIVITY, true, bind(&OptionsWindow::startPolygonizing, this)),
				mAlgo(algo)
			{
			}
			
			void startPolygonizing()
			{
				mAlgo.scheduleJob(bind(&MarchingCubes::polygonize, &mAlgo, mSlider.get()));
			}
		};

		void polygonize(size_t& value)
		{
			unique_lock<mutex> lock(mMutex);

			float isoValue = value * (maxIsoValue / SENSITIVITY);
			Color color = Color(1.0, 0.0, 0.0);

			if(debug == true)
			{
				for(int i = 0; i < platforms.size(); i++)
				{
					platforms[i].getInfo(CL_PLATFORM_NAME, &infoString);
					debugLog() << "platform #" << i << " name: " << infoString << endl;
					platforms[i].getInfo(CL_PLATFORM_VERSION, &infoString);
					debugLog() << "platform #" << i << " version: " << infoString << endl;
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
			}

			cl::Kernel kernel(program, "marchingCubes");

			// prepare buffers
			cl::Buffer bufferTriPoints = cl::Buffer(cont, CL_MEM_WRITE_ONLY, 12 * numCells * sizeof(cl_float4));

			cl::Buffer bufferIndices = cl::Buffer(cont, CL_MEM_WRITE_ONLY, numCells * sizeof(int));

			//cl::Buffer bufferFloatTest = cl::Buffer(context, CL_MEM_WRITE_ONLY, numCells * sizeof(float));
			//cl::Buffer bufferIntTest = cl::Buffer(context, CL_MEM_WRITE_ONLY, numCells * sizeof(int));


			kernel.setArg(0, isoValue);
			kernel.setArg(1, bufferEdgeTable);
			kernel.setArg(2, bufferTriTable);
			kernel.setArg(3, bufferValues);
			kernel.setArg(4, bufferPointsVec);
			kernel.setArg(5, bufferTriPoints);
			kernel.setArg(6, bufferIndices);

			//kernel.setArg(7, bufferFloatTest);
			//kernel.setArg(8, bufferIntTest);

			// run the kernel
			cl::NDRange global(numCells);
			cl::NDRange local(1);
			queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);

			// get the results
			/*
			float* floatTest = new float[numCells];
			queue.enqueueReadBuffer(bufferFloatTest, CL_TRUE, 0, numCells * sizeof(float), floatTest);
			for(int i = 0; i < numCells; i++)
			{
				debugLog() << floatTest[i] << endl;
			}

			int* intTest = new int[numCells];
			queue.enqueueReadBuffer(bufferIntTest, CL_TRUE, 0, numCells * sizeof(int), intTest);
			for(int i = 0; i < numCells; i++)
			{
				debugLog() << intTest[i] << endl;
			}
			*/
				
			cl_float4* triPoints = new cl_float4[12 * numCells];
			queue.enqueueReadBuffer(bufferTriPoints, CL_TRUE, 0, 12 * numCells * sizeof(cl_float4), triPoints);

			int* indices = new int[numCells];
			queue.enqueueReadBuffer(bufferIndices, CL_TRUE, 0, numCells * sizeof(int), indices);


			// draw the triangles
			polygonGroup = makeGraphics("iso surface");
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
			//delete[] floatTest;
			//delete[] intTest;
			delete[] triPoints;
			delete[] indices;
		}

		mutex mMutex;
		Window<OptionsWindow> mWindow;
		unique_ptr<Graphics> polygonGroup;

		bool debug;
		float maxIsoValue;

		size_t numCells;
		size_t numCellPoints;
		size_t numValues;
		float* values;
		cl_float4* pointsVec;

		cl_int error;
		cl_ulong infoNumber;
		string infoString;

		vector<cl::Platform> platforms;
		vector<cl::Device> devices;
		cl::Context cont;
		cl::CommandQueue queue;
		cl::Program program;
		int* triTableArray;

		cl::Buffer bufferEdgeTable;
		cl::Buffer bufferTriTable;
		cl::Buffer bufferValues;
		cl::Buffer bufferPointsVec;

		MarchingCubes(const Parameters& parameters): Algorithm(parameters), mWindow(*this)
		{
			// get the data
			const TensorFieldGridBased<3, Scalar>* field = parameters.get<const TensorFieldGridBased<3, Scalar>*>("field");
			const Grid<3>* grid = parameters.get<const Grid<3>*>("grid");
			maxIsoValue = parameters.get<float>("max iso value");
			debug = parameters.get<bool>("debug");
			

			if(field == false)
			{
				throw runtime_error("field not set!");
			}
			if(grid == false)
			{
				throw runtime_error("grid not set!");
			}
 
			auto discreteEvaluator = field->makeDiscreteEvaluator();
			auto& points = grid->parent().points();

			polygonGroup = makeGraphics("iso surface");

			numCells = grid->numCells();
			//numCells = 1000;
			numCellPoints = 8;
			numValues = numCells * numCellPoints;
			int time = 0;

			// load all scalar values into an array
			values = new float[numValues];
			for(Progress i(*this, "load values", numCells); i < numCells; ++i)
			{	
				Cell cell = grid->cell(i);

				values[(i * numCellPoints) + 0] = discreteEvaluator->value(cell.index(0))();
				values[(i * numCellPoints) + 1] = discreteEvaluator->value(cell.index(1))();
				values[(i * numCellPoints) + 2] = discreteEvaluator->value(cell.index(2))();
				values[(i * numCellPoints) + 3] = discreteEvaluator->value(cell.index(3))();

				values[(i * numCellPoints) + 4] = discreteEvaluator->value(cell.index(4))();
				values[(i * numCellPoints) + 5] = discreteEvaluator->value(cell.index(5))();
				values[(i * numCellPoints) + 6] = discreteEvaluator->value(cell.index(6))();
				values[(i * numCellPoints) + 7] = discreteEvaluator->value(cell.index(7))();
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

			triTableArray = new int[256 * 16];
			for(int i = 0; i < 256; ++i)
			{
				for(int j = 0; j < 16; ++j)
				{
					triTableArray[(i * 16) + j] = triTable[i][j];
				}
			}

			// prepare everything
			cl::Platform::get(&platforms);

			cl_context_properties properties[3] = { 
				CL_CONTEXT_PLATFORM, 
				(cl_context_properties)(platforms[0])(), 
				0 
			};
			cl::Context context(CL_DEVICE_TYPE_GPU, properties, NULL, NULL, &error);
			printError(error, "context");

			cont = context;

			devices = cont.getInfo<CL_CONTEXT_DEVICES>();
			if(devices.size() == 0)
			{
				throw runtime_error("No devices found!");
			}

			queue = cl::CommandQueue(context, devices[0]);

			// prepare the kernel
			string s = getKernelSource();
			cl::Program::Sources source(1, std::make_pair(s.c_str(), s.length()+1));

			program = cl::Program(cont, source);
			error = program.build(devices);
			printError(error, "build");

			program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &infoString);
			debugLog() << "build info: " << infoString << endl;

			// prepare buffers
			bufferEdgeTable = cl::Buffer(context, CL_MEM_READ_ONLY, 256 * sizeof(int));
			error = queue.enqueueWriteBuffer(bufferEdgeTable, CL_TRUE, 0, 256 * sizeof(int), edgeTable);
			printError(error, "bufferEdgeTable");

			bufferTriTable = cl::Buffer(context, CL_MEM_READ_ONLY, 256 * 16 * sizeof(int));
			error = queue.enqueueWriteBuffer(bufferTriTable, CL_TRUE, 0, 256 * 16 * sizeof(int), triTableArray);
			printError(error, "bufferTriTable");

			bufferValues = cl::Buffer(context, CL_MEM_READ_ONLY, numValues * sizeof(float));
			error = queue.enqueueWriteBuffer(bufferValues, CL_TRUE, 0, numValues * sizeof(float), values);
			printError(error, "bufferValues");

			bufferPointsVec = cl::Buffer(cont, CL_MEM_READ_ONLY, numValues *sizeof(cl_float4));
			error = queue.enqueueWriteBuffer(bufferPointsVec, CL_TRUE, 0, numValues *sizeof(cl_float4), pointsVec);
			printError(error, "bufferPointsVec");

			// polygonizes once at startup
			size_t startValue = 0;
			polygonize(startValue);
		}
	};

	AlgorithmRegister<MarchingCubes> dummy("MarchingCubesCL", "executes MarchingCubesCL");
}
